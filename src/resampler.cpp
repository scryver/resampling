#include "../libberdip/src/platform.h"
#include "../libberdip/src/std_memory.h"

#include "../xfloat/src/xfloat.h"
#include "../xfloat/src/xfloat_math.h"

#include "../sound-things/src/wav.h"

global API api;
global MemoryAPI *gMemoryApi = &api.memory;

#include "../libberdip/src/memory.cpp"
#include "../libberdip/src/std_memory.cpp"
#include "../libberdip/src/std_file.c"

#include "../libberdip/src/argoptions.cpp"

#include "./xfloat_constants_16.c"
#include "../xfloat/src/xfloat.c"
#include "../xfloat/src/xfloat_math.c"
#include "../xfloat/src/xfloat_string.c"

#include "../xfloat/src/float512.h"

#include "../sound-things/src/wav.cpp"

#include "coefgen.cpp"

struct Resampler
{
    u32 fsIn;
    u32 fsOut;
    u32 fsMain;
    
    u32 prediv;
    u32 L;
    u32 M;
    
    u32 coefCount;
    f64 *coefs;
    
    f64 coefMult;
    u32 coefOffset;
    u32 sampleStep;
    s32 sampleIdx;
    
    u32 delayCount;
    f64 *delayBuf;
};

internal Resampler *
create_resampler(MemoryArena *memory, u32 fin, u32 fout, u32 fmain, u32 coefCount, f64 *coefs)
{
    Resampler *result = arena_allocate_struct(memory, Resampler, default_memory_alloc());
    result->fsIn = fin;
    result->fsOut = fout;
    result->fsMain = fmain;
    
    result->coefCount = coefCount;
    result->coefs = coefs;
    
    result->prediv = 1;
    result->L = fmain / fin;
    result->M = fmain / fout;
    
    while (((result->L % 2) == 0) &&
           ((result->M % 2) == 0))
    {
        result->prediv *= 2;
        result->L /= 2;
        result->M /= 2;
    }
    
    if ((result->L < result->M) &&
        ((result->M % result->L) == 0))
    {
        // NOTE(michiel): Pure decimation / down sampling
        result->prediv *= result->L;
        result->M /= result->L;
        result->L = 1;
    }
    else if ((result->M < result->L) &&
             ((result->L % result->M) == 0))
    {
        // NOTE(michiel): Pure interpolation / up sampling
        result->prediv *= result->M;
        result->L /= result->M;
        result->M = 1;
    }
    else if (result->L == result->M)
    {
        fprintf(stderr, "No resampling needed!\n");
    }
    
    result->coefMult = (f64)result->prediv * (f64)result->L;
    result->coefOffset = 0;
    result->sampleStep = result->M / result->L;
    
    result->delayCount = coefCount / result->prediv + 1;
    result->delayBuf = arena_allocate_array(memory, f64, result->delayCount, default_memory_alloc());
    
    return result;
}

internal void
resample(Resampler *resampler, u32 inputCount, f64 *input, u32 outputCount, f64 *output)
{
    f64 coefMult = (f64)resampler->prediv * (f64)resampler->L;
    
    i_expect(outputCount >= ((inputCount * resampler->L) / resampler->M));
    
    u32 sampleStep = resampler->M / resampler->L;
    u32 coefOffset = 0;
    u32 sampleIdx  = 0;
    
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        u32 cIdx = coefOffset * resampler->prediv;
        s32 sIdx = sampleIdx;
        
        f64 sample = 0.0;
        while (cIdx < resampler->coefCount)
        {
            sample += coefMult * resampler->coefs[cIdx] * input[sIdx];
            cIdx += resampler->L * resampler->prediv;
            sIdx -= 1;
            if (sIdx < 0) {
                break;
            }
        }
        output[outIndex] = sample;
        
        u32 prevOffset = coefOffset;
        coefOffset = (coefOffset + resampler->M) % resampler->L;
        sampleIdx += sampleStep + ((prevOffset > coefOffset) ? 1 : 0);
    }
    i_expect(sampleIdx == inputCount);
}

internal void
resample_chunk(Resampler *resampler, u32 inputCount, f64 *input, u32 outputCount, f64 *output)
{
    i_expect(outputCount >= ((inputCount * resampler->L) / resampler->M));
    
    s32 sampleIdx  = resampler->sampleIdx;
    
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        u32 cIdx = resampler->coefOffset * resampler->prediv;
        s32 sIdx = sampleIdx;
        
        f64 sample = 0.0;
        while (cIdx < resampler->coefCount)
        {
            if (sIdx < 0) {
                u32 delayIdx = resampler->delayCount + sIdx;
                i_expect(delayIdx < resampler->delayCount);
                sample += resampler->coefMult * resampler->coefs[cIdx] * resampler->delayBuf[delayIdx];
            } else {
                sample += resampler->coefMult * resampler->coefs[cIdx] * input[sIdx];
            }
            cIdx += resampler->L * resampler->prediv;
            sIdx -= 1;
        }
        output[outIndex] = sample;
        
        u32 prevOffset = resampler->coefOffset;
        resampler->coefOffset = (resampler->coefOffset + resampler->M) % resampler->L;
        sampleIdx += resampler->sampleStep + ((prevOffset > resampler->coefOffset) ? 1 : 0);
    }
    
    resampler->sampleIdx = sampleIdx - inputCount;
    
    if (resampler->delayCount <= inputCount)
    {
        for (u32 index = 0; index < resampler->delayCount; ++index)
        {
            resampler->delayBuf[index] = input[inputCount - resampler->delayCount + index];
        }
    }
    else
    {
        u32 remain = resampler->delayCount - inputCount;
        for (u32 index = 0; index < remain; ++index)
        {
            resampler->delayBuf[index] = resampler->delayBuf[resampler->delayCount - remain + index];
        }
        for (u32 index = remain; index < resampler->delayCount; ++index)
        {
            resampler->delayBuf[index] = input[index - remain];
        }
    }
}

internal b32
resample_flush(Resampler *resampler, u32 outputCount, f64 *output)
{
    b32 completeFlush = false;
    u32 sampleIdx = resampler->sampleIdx;
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        u32 cIdx = resampler->coefOffset * resampler->prediv;
        s32 sIdx = sampleIdx;
        
        f64 sample = 0.0;
        while (cIdx < resampler->coefCount)
        {
            if (sIdx < 0) {
                u32 delayIdx = resampler->delayCount + sIdx;
                i_expect(delayIdx < resampler->delayCount);
                sample += resampler->coefMult * resampler->coefs[cIdx] * resampler->delayBuf[delayIdx];
            }
            cIdx += resampler->L * resampler->prediv;
            sIdx -= 1;
        }
        output[outIndex] = sample;
        
        u32 prevOffset = resampler->coefOffset;
        resampler->coefOffset = (resampler->coefOffset + resampler->M) % resampler->L;
        sampleIdx += resampler->sampleStep + ((prevOffset > resampler->coefOffset) ? 1 : 0);
    }
    
    resampler->sampleIdx = sampleIdx;
    
    if (resampler->delayCount <= sampleIdx)
    {
        for (u32 index = 0; index < resampler->delayCount; ++index)
        {
            resampler->delayBuf[index] = 0.0;
        }
        resampler->sampleIdx = 0;
        resampler->coefOffset = 0;
        completeFlush = true;
    }
    else
    {
        u32 remain = resampler->delayCount - sampleIdx;
        for (u32 index = 0; index < remain; ++index)
        {
            resampler->delayBuf[index] = resampler->delayBuf[resampler->delayCount - remain + index];
        }
        for (u32 index = remain; index < resampler->delayCount; ++index)
        {
            resampler->delayBuf[index] = 0.0;
        }
    }
    return completeFlush;
}

int main(int argc, char **argv)
{
    MemoryArena arena = {};
    std_memory_api(gMemoryApi);
    std_file_api(&api.file);
    
    ArgOption argOptions[] =
    {
        arg_option_bool('a', static_string("first"), false, false, static_string("First item.")),
        arg_option_bool('b', static_string("second"), true, false, static_string("Second item.")),
        arg_option_counter('c', static_string("third"), false, 0, static_string("Third item.")),
        arg_option_counter('d', static_string("fourth"), true, 1, static_string("Fourth item.")),
        arg_option_int('e', static_string("fifth"), false, 2, static_string("Fifth item.")),
        arg_option_int('f', static_string("sixth"), true,  3, static_string("Sixth item.")),
        arg_option_float('g', static_string("seventh"), false,  1.0, static_string("Seventh item.")),
        arg_option_float('h', static_string("eight"), true,  4.0, static_string("Eight item.")),
        arg_option_string('i', static_string("ninth"), false, static_string("None"), static_string("Ninth item.")),
        arg_option_string('j', static_string("tenth"), true, static_string("Some"), static_string("Tenth item.")),
        arg_option_unnamed(static_string("input"), false, static_string("0"), static_string("Eleventh item.")),
        arg_option_unnamed(static_string("output"), true, static_string("1"), static_string("Twelfth item.")),
    };
    
    if (parse_options(argc, argv, array_count(argOptions), argOptions))
    {
        fprintf(stderr, "Parse success\n");
        print_usage(argv[0], array_count(argOptions), argOptions);
    }
    else
    {
        fprintf(stderr, "Failed to parse options\n");
        print_usage(argv[0], array_count(argOptions), argOptions);
    }
    
    
    
    return 0;
    
    u32 coefCount = 500003;
    u32 sampleFreq = 56448000;
    u32 cutoffFreq = 20500;
    
    f64 *coefs64 = load_or_create_coefs(static_string("data/base"), coefCount, sampleFreq, cutoffFreq);
    
    u32 inputFreq  =  44100;
    u32 outputFreq = 192000;
    
    Resampler *resampler = create_resampler(&arena, inputFreq, outputFreq, sampleFreq, coefCount, coefs64);
    
    u32 inputCount = inputFreq * 10;
    f64 *input = arena_allocate_array(&arena, f64, inputCount, default_memory_alloc());
    
    u32 outputCount = (inputCount * resampler->L) / resampler->M;
    f64 *output = arena_allocate_array(&arena, f64, outputCount, default_memory_alloc());
    
    for (u32 idx = 0; idx < inputCount; ++idx)
    {
        f64 sample = (f64)idx * 2.0 * F64_PI * 220.0 / (f64)inputFreq;
        sample = 0.1 * sin_pi(sample);
        input[idx] = sample;
    }
    
    {
        u32 chunkSize = 4096;
        u32 chunkCount = inputCount / chunkSize;
        // NOTE(michiel): Ignore last bit for now
        u32 chunkOutSize = (chunkSize * resampler->L) / resampler->M;
        
        f64 *src = input;
        f64 *dst = output;
        for (u32 chunkIdx = 0; chunkIdx < chunkCount; ++chunkIdx)
        {
            i_expect(src < (input + inputCount));
            i_expect(dst < (output + outputCount));
            resample_chunk(resampler, chunkSize, src, chunkOutSize, dst);
            src += chunkSize;
            dst += chunkOutSize;
        }
        
        while (resample_flush(resampler, chunkOutSize, dst))
        {
            i_expect(dst < (output + outputCount));
        }
    }
    
    //resample(resampler, inputCount, input, outputCount, output);
    
    WavSettings inputWav = {};
    inputWav.channelCount = 2;
    inputWav.sampleFrequency = inputFreq;
    inputWav.sampleResolution = 16;
    inputWav.sampleFrameSize = 4;
    inputWav.format = WavFormat_PCM;
    
    WavSettings outputWav = {};
    outputWav.channelCount = 2;
    outputWav.sampleFrequency = outputFreq;
    outputWav.sampleResolution = 16;
    outputWav.sampleFrameSize = 4;
    outputWav.format = WavFormat_PCM;
    
    Buffer inputData;
    inputData.size = sizeof(s16) * 2 * inputCount;
    inputData.data = arena_allocate_array(&arena, u8, inputData.size, default_memory_alloc());
    
    Buffer outputData;
    outputData.size = sizeof(s16) * 2 * outputCount;
    outputData.data = arena_allocate_array(&arena, u8, outputData.size, default_memory_alloc());
    
    s16 *inputD = (s16 *)inputData.data;
    for (u32 index = 0; index < inputCount; ++index)
    {
        f64 sample = input[index];
        inputD[2 * index + 0] = inputD[2 * index + 1] = (s16)((f64)S16_MAX * sample);
    }
    
    s16 *outputD = (s16 *)outputData.data;
    for (u32 index = 0; index < outputCount; ++index)
    {
        f64 sample = output[index];
        outputD[2 * index + 0] = outputD[2 * index + 1] = (s16)((f64)S16_MAX * sample);
    }
    
    wav_write_file(static_string("input.wav"), &inputWav, inputData);
    wav_write_file(static_string("output.wav"), &outputWav, outputData);
    
    return 0;
}
