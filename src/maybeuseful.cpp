
internal void
sinc(u32 count, f512 *dest, f512 *scale)
{
    // NOTE(michiel): sinc(x) = sin(pi*x) / (pi*x)
    i_expect((count % 2) == 1);
    u32 halfCount = count / 2;
    
    xf_copy(F512_ELEMENT_COUNT, gXF_One, dest[halfCount].e);
    for (u32 index = 0; index < halfCount; ++index)
    {
        u32 tapA = halfCount + index + 1;
        u32 tapB = halfCount - index - 1;
        f512 *value = dest + tapA;
        f512 tap = F512(index + 1);
        xf_mul(F512_ELEMENT_COUNT, tap.e, gXF_Pi, tap.e);
        xf_mul(F512_ELEMENT_COUNT, tap.e, scale->e, tap.e);
        // tap = pi*tap
        
        xf_sin(F512_ELEMENT_COUNT, tap.e, value->e);
        // value = sin(scale * pi * tap)
        
        xf_div(F512_ELEMENT_COUNT, value->e, tap.e, value->e);
        // value = sin(scale * pi * tap) / (scale * pi * tap)
        
        dest[tapB] = *value;
    }
}


internal void
calculate_coefs(u32 coefCount, f512 *coefs, u32 sampleFreq, u32 cutoffFreq)
{
    f512 coefCountX = F512(coefCount);
    f512 gaussNum = F512(0.084);
    gaussNum *= coefCountX;
    
    f512 fs = F512(sampleFreq);
    f512 f = F512(cutoffFreq);
    
    f512 w0;
    xf_mul(F512_ELEMENT_COUNT, gXF_Two, gXF_Pi, w0.e);
    w0 *= f;
    w0 /= fs;
    
    f512 minHalf = F512(gXF_Half);
    minHalf = -minHalf;
    
    f512 sum = {};
#if 1
    u32 halfCoef = coefCount / 2;
    //xf_div(F512_ELEMENT_COUNT, w0.e, gXF_Pi, coefs[halfCoef].e);
    xf_copy(F512_ELEMENT_COUNT, gXF_One, coefs[halfCoef].e);
    sum += coefs[halfCoef];
    for (u32 index = 0; index < halfCoef; ++index)
    {
        u32 tapA = halfCoef + index + 1;
        u32 tapB = halfCoef - index - 1;
        
        f512 *value = coefs + tapA;
        f512 tap = F512(index + 1);
        
#if 0        
        f512 win = tap / gaussNum;
        xf_sqr(F512_ELEMENT_COUNT, win.e, win.e);
        xf_mul(F512_ELEMENT_COUNT, win.e, minHalf.e, win.e);
        xf_exp(F512_ELEMENT_COUNT, win.e, win.e);
        // win = exp(-0.5 * (tap/gaussNum)^2)
#endif
        
        xf_mul(F512_ELEMENT_COUNT, tap.e, w0.e, value->e);
        xf_sin(F512_ELEMENT_COUNT, value->e, value->e);
        // value = sin(tap*w0)
        xf_mul(F512_ELEMENT_COUNT, tap.e, gXF_Pi, tap.e);
        xf_div(F512_ELEMENT_COUNT, value->e, tap.e, value->e);
        // value = sin(tap*w0)/(tap*pi)
        
        //xf_mul(F512_ELEMENT_COUNT, win.e, value->e, value->e);
        coefs[tapB] = *value;
        
        sum += *value;
        sum += *value;
    }
#else
    for (u32 index = 0; index < coefCount; ++index)
    {
        f512 *value = coefs + index;
        s32 tapIndex = (s32)index - (s32)coefCount / 2;
        if (tapIndex == 0)
        {
            xf_div(F512_ELEMENT_COUNT, w0.e, gXF_Pi, value->e);
            sum += *value;
        }
        else
        {
            f512 tap = F512(tapIndex);
            
            f512 win = tap / gaussNum;
            xf_sqr(F512_ELEMENT_COUNT, win.e, win.e);
            xf_mul(F512_ELEMENT_COUNT, win.e, minHalf.e, win.e);
            xf_exp(F512_ELEMENT_COUNT, win.e, win.e);
            // win = exp(-0.5 * (tap/gaussNum)^2)
            
            xf_mul(F512_ELEMENT_COUNT, tap.e, w0.e, value->e);
            xf_sin(F512_ELEMENT_COUNT, value->e, value->e);
            // value = sin(tap*w0)
            xf_mul(F512_ELEMENT_COUNT, tap.e, gXF_Pi, tap.e);
            xf_div(F512_ELEMENT_COUNT, value->e, tap.e, value->e);
            // value = sin(tap*w0)/(tap*pi)
            
            xf_mul(F512_ELEMENT_COUNT, win.e, value->e, value->e);
            sum += *value;
        }
    }
#endif
    
#if 0    
    for (u32 index = 0; index < coefCount; ++index)
    {
        f512 *value = coefs + index;
        (*value) /= sum;
    }
#endif
    
}

{
    
    f512 *coefs = arena_allocate_array(&arena, f512, coefCount, default_memory_alloc());
    calculate_coefs(coefCount, coefs, sampleFreq, cutoffFreq);
    
#if 1
    fprintf(stdout, "[  1000] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[  1000].e, 12);
    fprintf(stdout, "\n[  8000] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[  8000].e, 12);
    fprintf(stdout, "\n[ 16000] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[ 16000].e, 12);
    fprintf(stdout, "\n[ 32000] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[ 32000].e, 12);
    fprintf(stdout, "\n[ 64000] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[ 64000].e, 12);
    fprintf(stdout, "\n[125000] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[125000].e, 12);
    fprintf(stdout, "\n[249999] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[249999].e, 12);
    fprintf(stdout, "\n[250000] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[250000].e, 12);
    fprintf(stdout, "\n[250001] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[250001].e, 12);
    fprintf(stdout, "\n[250002] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[250002].e, 12);
    fprintf(stdout, "\n[250003] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[250003].e, 12);
    fprintf(stdout, "\n[375002] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[375002].e, 12);
    fprintf(stdout, "\n[436002] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[436002].e, 12);
    fprintf(stdout, "\n[468002] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[468002].e, 12);
    fprintf(stdout, "\n[484002] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[484002].e, 12);
    fprintf(stdout, "\n[492002] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[492002].e, 12);
    fprintf(stdout, "\n[499002] = ");
    xf_print(F512_ELEMENT_COUNT, coefs[499002].e, 12);
    fprintf(stdout, "\n");
#endif
    
#if 1
    f64 *coefs64 = arena_allocate_array(&arena, f64, coefCount, default_memory_alloc());
    for (u32 index = 0; index < coefCount; ++index)
    {
        f512 *src = coefs + index;
        f64 *dst = coefs64 + index;
        *dst = f64_from_xf(F512_ELEMENT_COUNT, src->e);
    }
    
    fprintf(stdout, "[  1000] = %g\n", coefs64[  1000]);
    fprintf(stdout, "[  8000] = %g\n", coefs64[  8000]);
    fprintf(stdout, "[ 16000] = %g\n", coefs64[ 16000]);
    fprintf(stdout, "[ 32000] = %g\n", coefs64[ 32000]);
    fprintf(stdout, "[ 64000] = %g\n", coefs64[ 64000]);
    fprintf(stdout, "[125000] = %g\n", coefs64[125000]);
    fprintf(stdout, "[249999] = %g\n", coefs64[249999]);
    fprintf(stdout, "[250000] = %g\n", coefs64[250000]);
    fprintf(stdout, "[250001] = %g\n", coefs64[250001]);
    fprintf(stdout, "[250002] = %g\n", coefs64[250002]);
    fprintf(stdout, "[250003] = %g\n", coefs64[250003]);
    fprintf(stdout, "[375002] = %g\n", coefs64[375002]);
    fprintf(stdout, "[436002] = %g\n", coefs64[436002]);
    fprintf(stdout, "[468002] = %g\n", coefs64[468002]);
    fprintf(stdout, "[484002] = %g\n", coefs64[484002]);
    fprintf(stdout, "[492002] = %g\n", coefs64[492002]);
    fprintf(stdout, "[499002] = %g\n", coefs64[499002]);
    
#endif
    
}