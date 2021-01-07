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
    
    u32 channelCount;
    u32 delayCount;
    f64 *delayBuf;
};

struct AudioBuffer
{
    u32 channelCount;
    u32 sampleCount;
    f64 *samples; // NOTE(michiel): samples[sampleCount][channelCount]
};

internal Resampler *
create_resampler(MemoryArena *memory, u32 channelCount, u32 fin, u32 fout, u32 fmain, u32 coefCount, f64 *coefs)
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
    result->sampleIdx = 0;
    
    result->channelCount = channelCount;
    result->delayCount = channelCount * (coefCount / result->prediv + 1);
    result->delayBuf = arena_allocate_array(memory, f64, result->delayCount, default_memory_alloc());
    
    return result;
}

internal void
resample(Resampler *resampler, u32 inputCount, f64 *input, u32 outputCount, f64 *output)
{
    f64 coefMult = (f64)resampler->prediv * (f64)resampler->L;
    
    i_expect(outputCount >= ((inputCount * resampler->L) / resampler->M));
    
    u32 sampleStep = resampler->channelCount * resampler->sampleStep;
    u32 coefOffset = 0;
    u32 sampleIdx  = 0;
    
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        for (u32 chanIdx = 0; chanIdx < resampler->channelCount; ++chanIdx)
        {
            u32 cIdx = coefOffset * resampler->prediv;
            s32 sIdx = sampleIdx + chanIdx;
            
            f64 sample = 0.0;
            while (cIdx < resampler->coefCount)
            {
                sample += coefMult * resampler->coefs[cIdx] * input[sIdx];
                cIdx += resampler->L * resampler->prediv;
                sIdx -= resampler->channelCount;
                if (sIdx < 0) {
                    break;
                }
            }
            output[resampler->channelCount * outIndex + chanIdx] = sample;
        }
        
        u32 prevOffset = coefOffset;
        coefOffset = (coefOffset + resampler->M) % resampler->L;
        sampleIdx += sampleStep + ((prevOffset > coefOffset) ? resampler->channelCount : 0);
    }
    i_expect(sampleIdx == (inputCount * resampler->channelCount));
}

internal void
resample_chunk(Resampler *resampler, u32 inputCount, f64 *input, u32 outputCount, f64 *output)
{
    i_expect(outputCount >= ((inputCount * resampler->L) / resampler->M));
    
    u32 sampleStep = resampler->channelCount * resampler->sampleStep;
    s32 sampleIdx  = resampler->sampleIdx;
    
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        for (u32 chanIdx = 0; chanIdx < resampler->channelCount; ++chanIdx)
        {
            u32 cIdx = resampler->coefOffset * resampler->prediv;
            s32 sIdx = sampleIdx + chanIdx;
            
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
                sIdx -= resampler->channelCount;
            }
            output[resampler->channelCount * outIndex + chanIdx] = sample;
        }
        
        u32 prevOffset = resampler->coefOffset;
        resampler->coefOffset = (resampler->coefOffset + resampler->M) % resampler->L;
        sampleIdx += sampleStep + ((prevOffset > resampler->coefOffset) ? resampler->channelCount : 0);
    }
    
    u32 totalInputCount = inputCount * resampler->channelCount;
    resampler->sampleIdx = sampleIdx - totalInputCount;
    
    if (resampler->delayCount <= totalInputCount)
    {
        for (u32 index = 0; index < resampler->delayCount; ++index)
        {
            resampler->delayBuf[index] = input[totalInputCount - resampler->delayCount + index];
        }
    }
    else
    {
        u32 remain = resampler->delayCount - totalInputCount;
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
        for (u32 chanIdx = 0; chanIdx < resampler->channelCount; ++chanIdx)
        {
            u32 cIdx = resampler->coefOffset * resampler->prediv;
            s32 sIdx = sampleIdx + chanIdx;
            
            f64 sample = 0.0;
            while (cIdx < resampler->coefCount)
            {
                if (sIdx < 0) {
                    u32 delayIdx = resampler->delayCount + sIdx;
                    i_expect(delayIdx < resampler->delayCount);
                    sample += resampler->coefMult * resampler->coefs[cIdx] * resampler->delayBuf[delayIdx];
                }
                cIdx += resampler->L * resampler->prediv;
                sIdx -= resampler->channelCount;
            }
            output[resampler->channelCount * outIndex + chanIdx] = sample;
        }
        
        u32 prevOffset = resampler->coefOffset;
        resampler->coefOffset = (resampler->coefOffset + resampler->M) % resampler->L;
        sampleIdx += resampler->channelCount * (resampler->sampleStep + ((prevOffset > resampler->coefOffset) ? 1 : 0));
    }
    
    resampler->sampleIdx = sampleIdx;;
    
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

internal void
f64_from_s16(u32 count, s16 *input, f64 *output)
{
    f64 oneOverS16 = 1.0 / (f64)(S16_MAX);
    for (u32 index = 0; index < count; ++index)
    {
        output[index] = (f64)input[index] * oneOverS16;
    }
}

internal void
s16_from_f64(u32 count, f64 *input, s16 *output)
{
    f64 s16Gain = (f64)(S16_MAX);
    for (u32 index = 0; index < count; ++index)
    {
        output[index] = s16_from_f64_round(input[index] * s16Gain);
    }
}

internal void
f64_from_s24(u32 count, u8 *input, f64 *output)
{
    f64 oneOverS24 = 1.0 / (f64)((s32)0x7FFFFFFF);
    for (u32 index = 0; index < count; ++index)
    {
        s32 value = (s32)((input[3*index + 2] << 24) |
                          (input[3*index + 1] << 16) |
                          (input[3*index + 0] <<  8));
        output[index] = (f64)value * oneOverS24;
    }
}

internal void
s24_from_f64(u32 count, f64 *input, u8 *output)
{
    f64 s24Gain = (f64)((s32)0x007FFFFF);
    for (u32 index = 0; index < count; ++index)
    {
        s32 value = round64(input[index] * s24Gain);
        output[3*index + 0] = (value >>  0) & 0xFF;
        output[3*index + 1] = (value >>  8) & 0xFF;
        output[3*index + 2] = (value >> 16) & 0xFF;
    }
}

internal void
f64_from_s32(u32 count, s32 *input, f64 *output)
{
    f64 oneOverS32 = 1.0 / (f64)(S32_MAX);
    for (u32 index = 0; index < count; ++index)
    {
        output[index] = (f64)input[index] * oneOverS32;
    }
}

internal void
s32_from_f64(u32 count, f64 *input, s32 *output)
{
    f64 s32Gain = (f64)(S32_MAX);
    for (u32 index = 0; index < count; ++index)
    {
        output[index] = s32_from_f64_round(input[index] * s32Gain);
    }
}

int main(int argc, char **argv)
{
    MemoryArena arena = {};
    std_memory_api(gMemoryApi);
    std_file_api(&api.file);
    
#if 0
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
#endif
    
    u32 coefCount = 500003;
    u32 sampleFreq = 56448000;
    u32 cutoffFreq = 20500;
    
    f64 *coefs64 = load_or_create_coefs(static_string("data/base"), coefCount, sampleFreq, cutoffFreq);
    
    u32 channelCount = 2;
    u32 inputFreq  = 192000;
    u32 outputFreq =  44100;
    
    if (argc == 3)
    {
        String filename = string(argv[1]);
        String outFreqName = string(argv[2]);
        outputFreq = number_from_string(outFreqName);
        
        MemoryAllocator alloc = {};
        initialize_arena_allocator(&arena, &alloc);
        WavReader reader = wav_load_file(&alloc, filename);
        if (reader.settings.sampleFrequency)
        {
            inputFreq = reader.settings.sampleFrequency;
            channelCount = reader.settings.channelCount;
            
            Resampler *resampler = create_resampler(&arena, channelCount, inputFreq, outputFreq, sampleFreq, coefCount, coefs64);
            
            u32 inputCount = reader.dataCount / (reader.settings.sampleFrameSize);
            f64 *inputSamples = arena_allocate_array(&arena, f64, channelCount * inputCount, default_memory_alloc());
            
            u32 outputCount = (u32)(((u64)inputCount * (u64)resampler->L) / resampler->M);
            f64 *outputSamples = arena_allocate_array(&arena, f64, channelCount * outputCount, default_memory_alloc());
            
            switch (reader.settings.sampleResolution)
            {
                case 16: {
                    f64_from_s16(channelCount * inputCount, (s16 *)(reader.rawData.data + reader.dataOffset), inputSamples);
                } break;
                
                case 24: {
                    i_expect(reader.settings.sampleFrameSize == 3*channelCount);
                    f64_from_s24(channelCount * inputCount, reader.rawData.data + reader.dataOffset, inputSamples);
                } break;
                
                case 32: {
                    f64_from_s32(channelCount * inputCount, (s32 *)(reader.rawData.data + reader.dataOffset), inputSamples);
                } break;
                
                INVALID_DEFAULT_CASE;
            }
            
#if 0
            {
                u32 chunkSize = 4096;
                u32 chunkCount = inputCount / chunkSize;
                // NOTE(michiel): Ignore last bit for now
                u32 chunkOutSize = (chunkSize * resampler->L) / resampler->M;
                
                f64 *src = inputSamples;
                f64 *dst = outputSamples;
                for (u32 chunkIdx = 0; chunkIdx < chunkCount; ++chunkIdx)
                {
                    i_expect(src < (inputSamples + inputCount * resampler->channelCount));
                    i_expect(dst < (outputSamples + outputCount * resampler->channelCount));
                    resample_chunk(resampler, chunkSize, src, chunkOutSize, dst);
                    src += chunkSize * resampler->channelCount;
                    dst += chunkOutSize * resampler->channelCount;
                }
                
                while (!resample_flush(resampler, chunkOutSize, dst))
                {
                    i_expect(dst < (outputSamples + outputCount * resampler->channelCount));
                }
            }
#else
            resample(resampler, inputCount, inputSamples, outputCount, outputSamples);
#endif
            f64 maxSample = 0.0;
            f64 minSample = 0.0;
            
            s16 *outputData = arena_allocate_array(&arena, s16, channelCount * outputCount, default_memory_alloc());
            
            s16 *outputD = outputData;
            for (u32 index = 0; index < outputCount * channelCount; ++index)
            {
                f64 sample = outputSamples[index];
                outputD[index] = (s16)((f64)S16_MAX * sample);
                
                if (maxSample < sample) {
                    maxSample = sample;
                }
                if (minSample > sample) {
                    minSample = sample;
                }
            }
            fprintf(stderr, "Max: %f, min: %f\n", maxSample, minSample);
            
            maxSample = 0.0;
            minSample = 0.0;
            
            for (u32 index = 0; index < inputCount * channelCount; ++index)
            {
                f64 sample = inputSamples[index];
                if (maxSample < sample) {
                    maxSample = sample;
                }
                if (minSample > sample) {
                    minSample = sample;
                }
            }
            fprintf(stderr, "Max: %f, min: %f\n", maxSample, minSample);
            
            //s16_from_f64(channelCount * outputCount, outputSamples, outputData);
            
            Buffer outputBuffer;
            outputBuffer.size = sizeof(s16) * channelCount * outputCount;
            outputBuffer.data = (u8 *)outputData;
            
            WavSettings outputWav = {};
            outputWav.channelCount = channelCount;
            outputWav.sampleFrequency = outputFreq;
            outputWav.sampleResolution = 16;
            outputWav.sampleFrameSize = 4;
            outputWav.format = WavFormat_PCM;
            
            wav_write_file(static_string("resampled.wav"), &outputWav, outputBuffer);
            
#if 0            
            s16 *inputData = arena_allocate_array(&arena, s16, channelCount * inputCount, default_memory_alloc());
            s16_from_f64(channelCount * inputCount, inputSamples, inputData);
            
            Buffer inputBuffer;
            inputBuffer.size = sizeof(s16) * channelCount * inputCount;
            inputBuffer.data = (u8 *)inputData;
            
            WavSettings inputWav = {};
            inputWav.channelCount = channelCount;
            inputWav.sampleFrequency = inputFreq;
            inputWav.sampleResolution = 16;
            inputWav.sampleFrameSize = 4;
            inputWav.format = WavFormat_PCM;
            
            wav_write_file(static_string("presampled.wav"), &inputWav, inputBuffer);
#endif
        }
    }
    else
    {
        Resampler *resampler = create_resampler(&arena, channelCount, inputFreq, outputFreq, sampleFreq, coefCount, coefs64);
        
        u32 inputCount = inputFreq * 10;
        f64 *input = arena_allocate_array(&arena, f64, channelCount * inputCount, default_memory_alloc());
        
        u32 outputCount = (inputCount * resampler->L) / resampler->M;
        f64 *output = arena_allocate_array(&arena, f64, channelCount * outputCount, default_memory_alloc());
        
        for (u32 idx = 0; idx < inputCount; ++idx)
        {
            for (u32 chanIdx = 0; chanIdx < channelCount; ++chanIdx)
            {
                f64 sample = (f64)idx * 2.0 * F64_PI * (f64)(chanIdx + 1) * 220.0 / (f64)inputFreq;
                sample = 0.1 * sin_pi(sample);
                input[channelCount * idx + chanIdx] = sample;
            }
        }
        
#if 1
        {
            u32 chunkSize = 4096;
            u32 chunkCount = inputCount / chunkSize;
            // NOTE(michiel): Ignore last bit for now
            u32 chunkOutSize = (chunkSize * resampler->L) / resampler->M;
            
            f64 *src = input;
            f64 *dst = output;
            for (u32 chunkIdx = 0; chunkIdx < chunkCount; ++chunkIdx)
            {
                i_expect(src < (input + inputCount * resampler->channelCount));
                i_expect(dst < (output + outputCount * resampler->channelCount));
                resample_chunk(resampler, chunkSize, src, chunkOutSize, dst);
                src += chunkSize * resampler->channelCount;
                dst += chunkOutSize * resampler->channelCount;
            }
            
            while (resample_flush(resampler, chunkOutSize, dst))
            {
                i_expect(dst < (output + outputCount));
            }
        }
#else
        resample(resampler, inputCount, input, outputCount, output);
#endif
        
        WavSettings inputWav = {};
        inputWav.channelCount = channelCount;
        inputWav.sampleFrequency = inputFreq;
        inputWav.sampleResolution = 16;
        inputWav.sampleFrameSize = 4;
        inputWav.format = WavFormat_PCM;
        
        WavSettings outputWav = {};
        outputWav.channelCount = channelCount;
        outputWav.sampleFrequency = outputFreq;
        outputWav.sampleResolution = 16;
        outputWav.sampleFrameSize = 4;
        outputWav.format = WavFormat_PCM;
        
        Buffer inputData;
        inputData.size = sizeof(s16) * channelCount * inputCount;
        inputData.data = arena_allocate_array(&arena, u8, inputData.size, default_memory_alloc());
        
        Buffer outputData;
        outputData.size = sizeof(s16) * channelCount * outputCount;
        outputData.data = arena_allocate_array(&arena, u8, outputData.size, default_memory_alloc());
        
        s16 *inputD = (s16 *)inputData.data;
        for (u32 index = 0; index < inputCount * channelCount; ++index)
        {
            f64 sample = input[index];
            inputD[index] = (s16)((f64)S16_MAX * sample);
        }
        
        f64 maxSample = 0.0;
        f64 minSample = 0.0;
        
        s16 *outputD = (s16 *)outputData.data;
        for (u32 index = 0; index < outputCount * channelCount; ++index)
        {
            f64 sample = output[index];
            outputD[index] = (s16)((f64)S16_MAX * sample);
            
            if (maxSample < sample) {
                maxSample = sample;
            }
            if (minSample > sample) {
                minSample = sample;
            }
        }
        fprintf(stderr, "Max: %f, min: %f\n", maxSample, minSample);
        
        wav_write_file(static_string("input.wav"), &inputWav, inputData);
        wav_write_file(static_string("output.wav"), &outputWav, outputData);
    }
    
    return 0;
}
