#define DO_CHUNKED         0
#define DO_STEREO_OPT      1

#define CHUNK_SIZE         8192

#include <time.h>

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

internal struct timespec
linux_get_wall_clock(void)
{
    struct timespec clock;
    clock_gettime(CLOCK_MONOTONIC, &clock);
    return clock;
}

internal f32
linux_get_seconds_elapsed(struct timespec start, struct timespec end)
{
    return ((f32)(end.tv_sec - start.tv_sec)
            + ((f32)(end.tv_nsec - start.tv_nsec) * 1e-9f));
}

struct Resampler
{
    u32 fsIn;
    u32 fsOut;
    u32 fsMain;
    
    u32 prediv;
    u32 L;
    u32 M;
    
    u32 coefCount;
    f64 **coefs; // NOTE(michiel): coefs[L][coefCount]
    
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
        result->prediv *= result->M;
        result->L = 1;
        result->M = 1;
    }
    
    result->coefCount = coefCount / (result->L * result->prediv) + 1;
    i_expect(result->coefCount > 100);
    
    result->coefs = arena_allocate_array(memory, f64 *, result->L, default_memory_alloc());
    for (u32 tableIdx = 0; tableIdx < result->L; ++tableIdx)
    {
        result->coefs[tableIdx] = arena_allocate_array(memory, f64, result->coefCount, default_memory_alloc());
    }
    
    f64 coefMult = (f64)result->prediv * (f64)result->L;
    for (u32 index = 0; index < coefCount; index += result->prediv)
    {
        u32 divIdx = index / result->prediv;
        u32 tableIdx = divIdx % result->L;
        u32 coefIdx = divIdx / result->L;
        i_expect(coefIdx < result->coefCount);
        result->coefs[tableIdx][coefIdx] = coefMult * coefs[index];
    }
    
    for (u32 tableIdx = 0; tableIdx < result->L; ++tableIdx)
    {
        f64 summing = 0.0;
        for (u32 idx = 0; idx < result->coefCount; ++idx)
        {
            summing += result->coefs[tableIdx][idx];
        }
        //fprintf(stdout, "[%3u]: %f\n", tableIdx, summing);
        
        for (u32 idx = 0; idx < result->coefCount; ++idx)
        {
            result->coefs[tableIdx][idx] /= summing;
        }
    }
    
    result->coefOffset = 0;
    result->sampleStep = result->M / result->L;
    result->sampleIdx = 0;
    
    result->channelCount = channelCount;
    result->delayCount = channelCount * (result->L * result->coefCount);
    result->delayBuf = arena_allocate_array(memory, f64, result->delayCount, default_memory_alloc());
    
    return result;
}

internal void
resample(Resampler *resampler, u32 inputCount, f64 *input, u32 outputCount, f64 *output)
{
    i_expect(outputCount >= ((inputCount * resampler->L) / resampler->M));
    
    u32 sampleStep = resampler->channelCount * resampler->sampleStep;
    u32 coefOffset = 0;
    u32 sampleIdx  = 0;
    
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        f64 *coefs = resampler->coefs[coefOffset];
        
        for (u32 chanIdx = 0; chanIdx < resampler->channelCount; ++chanIdx)
        {
            s32 sIdx = sampleIdx + chanIdx;
            
            f64 sample = 0.0;
            for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
            {
                sample += coefs[cIdx] * input[sIdx];
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
resample_2ch_interleaved(Resampler *resampler, u32 inputCount, f64 *input, u32 outputCount, f64 *output)
{
    // NOTE(michiel): Idee om te proberen is de coefOffset inverteren en dan beide arrays laten optellen
    i_expect(resampler->channelCount == 2);
    i_expect(outputCount >= ((inputCount * resampler->L) / resampler->M));
    i_expect(((umm)input & 0xF) == 0);
    i_expect(((umm)output & 0xF) == 0);
    
    u32 sampleStep = 2 * resampler->sampleStep;
    u32 coefOffset = 0;
    u32 sampleIdx  = 0;
    
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        f64 *coefs = resampler->coefs[coefOffset];
        s32 sIdx = sampleIdx;
        
#if 0
        if (sIdx < 2 * resampler->coefCount)
        {
#if 1
            __m128d sample = _mm_set1_pd(0.0);
            f64 *src = input + sIdx;
            for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
            {
                __m128d coef = _mm_set1_pd(coefs[cIdx]);
                __m128d inp  = _mm_load_pd(src);
                sample = _mm_add_pd(sample, _mm_mul_pd(coef, inp));
                src -= 2;
                if (src < input) {
                    break;
                }
            }
            _mm_store_pd(output + 2 * outIndex, sample);
#else
            f64 sample0 = 0.0;
            f64 sample1 = 0.0;
            for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
            {
                sample0 += coefs[cIdx] * input[sIdx+0];
                sample1 += coefs[cIdx] * input[sIdx+1];
                sIdx -= 2;
                if (sIdx < 0) {
                    break;
                }
            }
            output[2 * outIndex + 0] = sample0;
            output[2 * outIndex + 1] = sample1;
#endif
        }
        else
        {
#if 1
            __m128d sample = _mm_set1_pd(0.0);
            f64 *src = input + sIdx;
            for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
            {
                __m128d coef = _mm_set1_pd(coefs[cIdx]);
                __m128d inp  = _mm_load_pd(src);
                sample = _mm_add_pd(sample, _mm_mul_pd(coef, inp));
                src -= 2;
            }
            _mm_store_pd(output + 2 * outIndex, sample);
#else
            f64 sample0 = 0.0;
            f64 sample1 = 0.0;
            for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
            {
                sample0 += coefs[cIdx] * input[sIdx+0];
                sample1 += coefs[cIdx] * input[sIdx+1];
                sIdx -= 2;
            }
            output[2 * outIndex + 0] = sample0;
            output[2 * outIndex + 1] = sample1;
#endif
        }
#else
        if (sIdx < 2 * resampler->coefCount)
        {
            f64 sample0 = 0.0;
            f64 sample1 = 0.0;
            for (u32 cIdx = 0; cIdx < resampler->coefCount / 4; ++cIdx)
            {
                sample0 += coefs[4*cIdx + 0] * input[sIdx+0];
                sample1 += coefs[4*cIdx + 0] * input[sIdx+1];
                sIdx -= 2;
                if (sIdx < 0) {
                    break;
                }
                sample0 += coefs[4*cIdx + 1] * input[sIdx+0];
                sample1 += coefs[4*cIdx + 1] * input[sIdx+1];
                sIdx -= 2;
                if (sIdx < 0) {
                    break;
                }
                sample0 += coefs[4*cIdx + 2] * input[sIdx+0];
                sample1 += coefs[4*cIdx + 2] * input[sIdx+1];
                sIdx -= 2;
                if (sIdx < 0) {
                    break;
                }
                sample0 += coefs[4*cIdx + 3] * input[sIdx+0];
                sample1 += coefs[4*cIdx + 3] * input[sIdx+1];
                sIdx -= 2;
                if (sIdx < 0) {
                    break;
                }
            }
            for (u32 cIdx = 0; cIdx < resampler->coefOffset % 4; ++cIdx)
            {
                sample0 += coefs[resampler->coefCount - ((resampler->coefOffset % 4) - cIdx)] * input[sIdx+0];
                sample1 += coefs[resampler->coefCount - ((resampler->coefOffset % 4) - cIdx)] * input[sIdx+1];
                sIdx -= 2;
                if (sIdx < 0) {
                    break;
                }
            }
            output[2 * outIndex + 0] = sample0;
            output[2 * outIndex + 1] = sample1;
        }
        else
        {
            __m128d sample = _mm_set1_pd(0.0);
            f64 *src = input + sIdx;
            for (u32 cIdx = 0; cIdx < resampler->coefCount / 4; ++cIdx)
            {
                __m128d coef = _mm_set1_pd(coefs[4*cIdx + 0]);
                __m128d inp  = _mm_load_pd(src + 0);
                sample = _mm_add_pd(sample, _mm_mul_pd(coef, inp));
                coef = _mm_set1_pd(coefs[4*cIdx + 1]);
                inp  = _mm_load_pd(src - 2);
                sample = _mm_add_pd(sample, _mm_mul_pd(coef, inp));
                coef = _mm_set1_pd(coefs[4*cIdx + 2]);
                inp  = _mm_load_pd(src - 4);
                sample = _mm_add_pd(sample, _mm_mul_pd(coef, inp));
                coef = _mm_set1_pd(coefs[4*cIdx + 3]);
                inp  = _mm_load_pd(src - 6);
                sample = _mm_add_pd(sample, _mm_mul_pd(coef, inp));
                src -= 8;
            }
            for (u32 cIdx = 0; cIdx < resampler->coefOffset % 4; ++cIdx)
            {
                __m128d coef = _mm_set1_pd(coefs[resampler->coefCount - ((resampler->coefOffset % 4) - cIdx)]);
                __m128d inp  = _mm_load_pd(src);
                sample = _mm_add_pd(sample, _mm_mul_pd(coef, inp));
                src -= 2;
            }
            _mm_store_pd(output + 2 * outIndex, sample);
            
#if 0            
            f64 sample0 = 0.0;
            f64 sample1 = 0.0;
            for (u32 cIdx = 0; cIdx < resampler->coefCount / 4; ++cIdx)
            {
                sample0 += coefs[4*cIdx + 0] * input[sIdx+0];
                sample1 += coefs[4*cIdx + 0] * input[sIdx+1];
                sample0 += coefs[4*cIdx + 1] * input[sIdx-2];
                sample1 += coefs[4*cIdx + 1] * input[sIdx-1];
                sample0 += coefs[4*cIdx + 2] * input[sIdx-4];
                sample1 += coefs[4*cIdx + 2] * input[sIdx-3];
                sample0 += coefs[4*cIdx + 3] * input[sIdx-6];
                sample1 += coefs[4*cIdx + 3] * input[sIdx-5];
                sIdx -= 8;
            }
            for (u32 cIdx = 0; cIdx < resampler->coefOffset % 4; ++cIdx)
            {
                sample0 += coefs[resampler->coefCount - ((resampler->coefOffset % 4) - cIdx)] * input[sIdx+0];
                sample1 += coefs[resampler->coefCount - ((resampler->coefOffset % 4) - cIdx)] * input[sIdx+1];
                sIdx -= 2;
            }
            output[2 * outIndex + 0] = sample0;
            output[2 * outIndex + 1] = sample1;
#endif
            
        }
#endif
        
        u32 prevOffset = coefOffset;
        coefOffset = (coefOffset + resampler->M) % resampler->L;
        sampleIdx += sampleStep + ((prevOffset > coefOffset) ? 2 : 0);
    }
    i_expect(sampleIdx == (inputCount * 2));
}

internal void
resample_chunk(Resampler *resampler, u32 inputCount, f64 *input, u32 outputCount, f64 *output)
{
    i_expect(outputCount >= ((inputCount * resampler->L) / resampler->M));
    
    u32 sampleStep = resampler->channelCount * resampler->sampleStep;
    s32 sampleIdx  = resampler->sampleIdx;
    
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        f64 *coefs = resampler->coefs[resampler->coefOffset];
        
        for (u32 chanIdx = 0; chanIdx < resampler->channelCount; ++chanIdx)
        {
            s32 sIdx = sampleIdx + chanIdx;
            
            f64 sample = 0.0;
            for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
            {
                if (sIdx < 0) {
                    u32 delayIdx = resampler->delayCount + sIdx;
                    i_expect(delayIdx < resampler->delayCount);
                    sample += coefs[cIdx] * resampler->delayBuf[delayIdx];
                } else {
                    sample += coefs[cIdx] * input[sIdx];
                }
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
        f64 *coefs = resampler->coefs[resampler->coefOffset];
        
        for (u32 chanIdx = 0; chanIdx < resampler->channelCount; ++chanIdx)
        {
            s32 sIdx = sampleIdx + chanIdx;
            
            f64 sample = 0.0;
            for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
            {
                if (sIdx < 0) {
                    u32 delayIdx = resampler->delayCount + sIdx;
                    i_expect(delayIdx < resampler->delayCount);
                    sample += coefs[cIdx] * resampler->delayBuf[delayIdx];
                }
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
resample_chunk_2ch_interleaved(Resampler *resampler, u32 inputCount, f64 *input, u32 outputCount, f64 *output)
{
    i_expect(resampler->channelCount == 2);
    i_expect(outputCount >= ((inputCount * resampler->L) / resampler->M));
    
    u32 sampleStep = 2 * resampler->sampleStep;
    s32 sampleIdx  = resampler->sampleIdx;
    
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        f64 *coefs = resampler->coefs[resampler->coefOffset];
        s32 sIdx = sampleIdx;
        
        f64 sample0 = 0.0;
        f64 sample1 = 0.0;
        for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
        {
            if (sIdx < 0) {
                i_expect(sIdx < -1);
                u32 delayIdx = resampler->delayCount + sIdx;
                i_expect(delayIdx < resampler->delayCount);
                sample0 += coefs[cIdx] * resampler->delayBuf[delayIdx + 0];
                sample1 += coefs[cIdx] * resampler->delayBuf[delayIdx + 1];
            } else {
                sample0 += coefs[cIdx] * input[sIdx + 0];
                sample1 += coefs[cIdx] * input[sIdx + 1];
            }
            sIdx -= 2;
            
            output[2 * outIndex + 0] = sample0;
            output[2 * outIndex + 1] = sample1;
        }
        
        u32 prevOffset = resampler->coefOffset;
        resampler->coefOffset = (resampler->coefOffset + resampler->M) % resampler->L;
        sampleIdx += sampleStep + ((prevOffset > resampler->coefOffset) ? 2 : 0);
    }
    
    u32 totalInputCount = inputCount * 2;
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
resample_flush_2ch_interleaved(Resampler *resampler, u32 outputCount, f64 *output)
{
    i_expect(resampler->channelCount == 2);
    b32 completeFlush = false;
    u32 sampleIdx = resampler->sampleIdx;
    for (u32 outIndex = 0; outIndex < outputCount; ++outIndex)
    {
        f64 *coefs = resampler->coefs[resampler->coefOffset];
        s32 sIdx = sampleIdx;
        
        f64 sample0 = 0.0;
        f64 sample1 = 0.0;
        for (u32 cIdx = 0; cIdx < resampler->coefCount; ++cIdx)
        {
            if (sIdx < 0) {
                i_expect(sIdx < -1);
                u32 delayIdx = resampler->delayCount + sIdx;
                i_expect(delayIdx < resampler->delayCount);
                sample0 += coefs[cIdx] * resampler->delayBuf[delayIdx + 0];
                sample1 += coefs[cIdx] * resampler->delayBuf[delayIdx + 1];
            }
            sIdx -= resampler->channelCount;
        }
        output[2 * outIndex + 0] = sample0;
        output[2 * outIndex + 1] = sample1;
        
        u32 prevOffset = resampler->coefOffset;
        resampler->coefOffset = (resampler->coefOffset + resampler->M) % resampler->L;
        sampleIdx += 2 * (resampler->sampleStep + ((prevOffset > resampler->coefOffset) ? 1 : 0));
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
    
#if 0    
    u32 freqList[] = {44100, 48000, 88200, 96000, 176400, 192000, 352800, 384000};
    for (u32 inIdx = 0; inIdx < array_count(freqList); ++inIdx)
    {
        for (u32 outIdx = 0; outIdx < array_count(freqList); ++outIdx)
        {
            u32 inFreq = freqList[inIdx];
            u32 outFreq = freqList[outIdx];
            
            Resampler *resampler = create_resampler(&arena, 2, inFreq, outFreq, sampleFreq, coefCount, coefs64);
            
            fprintf(stdout, "%6u -> %6u | coef: %4u x %4u\n", inFreq, outFreq, resampler->L, resampler->coefCount);
        }
    }
    return 0;
#endif
    
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
            
#if DO_CHUNKED
            {
                u32 chunkSize = CHUNK_SIZE;
                u32 chunkCount = inputCount / chunkSize;
                // NOTE(michiel): Ignore last bit for now
                u32 chunkOutSize = (chunkSize * resampler->L) / resampler->M;
                
                f64 *src = inputSamples;
                f64 *dst = outputSamples;
#if DO_STEREO_OPT
                if (resampler->channelCount == 2)
                {
                    for (u32 chunkIdx = 0; chunkIdx < chunkCount; ++chunkIdx)
                    {
                        i_expect(src < (inputSamples + inputCount * resampler->channelCount));
                        i_expect(dst < (outputSamples + outputCount * resampler->channelCount));
                        resample_chunk_2ch_interleaved(resampler, chunkSize, src, chunkOutSize, dst);
                        src += chunkSize * resampler->channelCount;
                        dst += chunkOutSize * resampler->channelCount;
                    }
                    
                    while (!resample_flush_2ch_interleaved(resampler, chunkOutSize, dst))
                    {
                        i_expect(dst < (outputSamples + outputCount * resampler->channelCount));
                    }
                }
                else
                {
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
#endif
            }
#else
            
#if DO_STEREO_OPT
            if (resampler->channelCount == 2) {
                resample_2ch_interleaved(resampler, inputCount, inputSamples, outputCount, outputSamples);
            } else {
                resample(resampler, inputCount, inputSamples, outputCount, outputSamples);
            }
#else
            resample(resampler, inputCount, inputSamples, outputCount, outputSamples);
#endif
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
    else if (argc == 2)
    {
        fprintf(stdout, "Doing tests\n");
        u32 testSeconds = 2;
        channelCount = 2;
        u32 freqList[] = {44100, 48000, 88200, 96000, 176400, 192000, 352800, 384000};
        f64 totalSecs = 0.0;
        u64 totalBytes = 0;
        for (u32 inIdx = 0; inIdx < array_count(freqList); ++inIdx)
        {
            for (u32 outIdx = 0; outIdx < array_count(freqList); ++outIdx)
            {
                TempArenaMemory tempMem = begin_temporary_memory(&arena);
                
                inputFreq = freqList[inIdx];
                outputFreq = freqList[outIdx];
                
                Resampler *resampler = create_resampler(&arena, channelCount, inputFreq, outputFreq, sampleFreq, coefCount, coefs64);
                
                u32 inputCount = inputFreq * testSeconds;
                f64 *inputSamples = arena_allocate_array(&arena, f64, channelCount * inputCount, align_memory_alloc(16));
                
                u32 outputCount = (inputCount * resampler->L) / resampler->M;
                f64 *outputSamples = arena_allocate_array(&arena, f64, channelCount * outputCount, align_memory_alloc(16));
                
                for (u32 idx = 0; idx < inputCount; ++idx)
                {
                    for (u32 chanIdx = 0; chanIdx < channelCount; ++chanIdx)
                    {
                        f64 sample = (f64)idx * 2.0 * F64_PI * (f64)(chanIdx + 1) * 220.0 / (f64)inputFreq;
                        sample = 0.1 * sin_pi(sample);
                        inputSamples[channelCount * idx + chanIdx] = sample;
                    }
                }
                
                fprintf(stdout, "%6u -> %6u | coef: %4u x %4u | ", inputFreq, outputFreq, resampler->L, resampler->coefCount);
                
                struct timespec startTime = linux_get_wall_clock();
#if DO_CHUNKED
                {
                    u32 chunkSize = CHUNK_SIZE;
                    u32 chunkCount = inputCount / chunkSize;
                    // NOTE(michiel): Ignore last bit for now
                    u32 chunkOutSize = (chunkSize * resampler->L) / resampler->M;
                    
                    f64 *src = inputSamples;
                    f64 *dst = outputSamples;
#if DO_STEREO_OPT
                    if (resampler->channelCount == 2)
                    {
                        for (u32 chunkIdx = 0; chunkIdx < chunkCount; ++chunkIdx)
                        {
                            i_expect(src < (inputSamples + inputCount * resampler->channelCount));
                            i_expect(dst < (outputSamples + outputCount * resampler->channelCount));
                            resample_chunk_2ch_interleaved(resampler, chunkSize, src, chunkOutSize, dst);
                            src += chunkSize * resampler->channelCount;
                            dst += chunkOutSize * resampler->channelCount;
                        }
                        
                        while (!resample_flush_2ch_interleaved(resampler, chunkOutSize, dst))
                        {
                            i_expect(dst < (outputSamples + outputCount * resampler->channelCount));
                        }
                    }
                    else
                    {
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
#endif
                }
#else
                
#if DO_STEREO_OPT
                if (resampler->channelCount == 2) {
                    resample_2ch_interleaved(resampler, inputCount, inputSamples, outputCount, outputSamples);
                } else {
                    resample(resampler, inputCount, inputSamples, outputCount, outputSamples);
                }
#else
                resample(resampler, inputCount, inputSamples, outputCount, outputSamples);
#endif
#endif
                f32 timing = linux_get_seconds_elapsed(startTime, linux_get_wall_clock());
                u64 bytes = (u64)resampler->L * resampler->coefCount * sizeof(f64);
                fprintf(stdout, "%6u: %f sec (%f usec/sample) (%lu kbytes)\n", outputFreq, timing, 1000000.0 * timing / (f64)outputCount, bytes / 1024);
                totalBytes += bytes;
                totalSecs += timing;
                
                end_temporary_memory(tempMem);
            }
        }
        
        fprintf(stdout, "\nTotal time: %f\n", totalSecs);
        char *suffix = "B";
        f64 totalBytesF = (f64)totalBytes;
        if (totalBytesF > 1024.0) {
            suffix = "kB";
            totalBytesF /= 1024.0;
        }
        if (totalBytesF > 1024.0) {
            suffix = "MB";
            totalBytesF /= 1024.0;
        }
        if (totalBytesF > 1024.0) {
            suffix = "GB";
            totalBytesF /= 1024.0;
        }
        fprintf(stdout, "Total bytes: %7.3f%s\n", totalBytesF, suffix);
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
        
#if DO_CHUNKED
        {
            u32 chunkSize = CHUNK_SIZE;
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
