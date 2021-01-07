struct CoefFile
{
    u32 fileMagic;
    u32 coefCount;
};

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
    
    //sinc(coefCount, coefs, &w0);
    
    f512 minHalf = F512(gXF_Half);
    minHalf = -minHalf;
    
    u32 halfCoef = coefCount / 2;
    
    f512 sum = {};
    sum += coefs[halfCoef];
    
    for (u32 index = 0; index < halfCoef; ++index)
    {
        u32 tapA = halfCoef + index + 1;
        u32 tapB = halfCoef - index - 1;
        
        f512 *value = coefs + tapA;
        f512 tap = F512(index + 1);
        
        f512 win = tap / gaussNum;
        xf_sqr(F512_ELEMENT_COUNT, win.e, win.e);
        xf_mul(F512_ELEMENT_COUNT, win.e, minHalf.e, win.e);
        xf_exp(F512_ELEMENT_COUNT, win.e, win.e);
        // win = exp(-0.5 * (tap/gaussNum)^2)
        
        xf_mul(F512_ELEMENT_COUNT, tap.e, w0.e, tap.e);
        // tap = scale * pi * tap
        
        xf_sin(F512_ELEMENT_COUNT, tap.e, value->e);
        // value = sin(scale * pi * tap)
        
        xf_div(F512_ELEMENT_COUNT, value->e, tap.e, value->e);
        // value = sin(scale * pi * tap) / (scale * pi * tap)
        
        xf_mul(F512_ELEMENT_COUNT, win.e, value->e, value->e);
        coefs[tapB] = *value;
        
        sum += *value;
        sum += *value;
    }
    
    for (u32 index = 0; index < coefCount; ++index)
    {
        f512 *value = coefs + index;
        (*value) /= sum;
    }
}

internal f64 *
load_or_create_coefs(String basename, u32 coefCount, u32 sampleFreq, u32 cutoffFreq)
{
    f64 *result = 0;
    
    u8 nameBuf[4096];
    String name = string_fmt(array_count(nameBuf), nameBuf, "%.*s_%u_fs%uHz_fc%uHz.cfs", STR_FMT(basename), coefCount, sampleFreq, cutoffFreq);
    
    MemoryAllocator alloc = {};
    initialize_platform_allocator(0, &alloc);
    Buffer readFile = api.file.read_entire_file(&alloc, name);
    if (readFile.size)
    {
        CoefFile *fileHdr = (CoefFile *)readFile.data;
        if ((fileHdr->fileMagic == MAKE_MAGIC('c', 'o', 'e', 'f')) &&
            (fileHdr->coefCount == coefCount))
        {
            result = (f64 *)(fileHdr + 1);
        }
        else
        {
            fprintf(stderr, "Invalid file\n");
        }
    }
    else
    {
        f512 *coefs = allocate_array(&alloc, f512, coefCount, default_memory_alloc());
        calculate_coefs(coefCount, coefs, sampleFreq, cutoffFreq);
        
        f64 *coefs64 = allocate_array(&alloc, f64, coefCount, default_memory_alloc());
        for (u32 index = 0; index < coefCount; ++index)
        {
            f512 *src = coefs + index;
            f64 *dst = coefs64 + index;
            *dst = f64_from_xf(F512_ELEMENT_COUNT, src->e);
        }
        
        CoefFile header = {};
        header.fileMagic = MAKE_MAGIC('c', 'o', 'e', 'f');
        header.coefCount = coefCount;
        
        ApiFile writeBack = api.file.open_file(name, FileOpen_Write);
        api.file.write_to_file(&writeBack, sizeof(CoefFile), &header);
        api.file.write_to_file(&writeBack, coefCount * sizeof(f64), coefs64);
        
        result = coefs64;
    }
    
    f64 sum = 0.0;
    for (u32 i = 0; i < coefCount; ++i)
    {
        sum += result[i];
    }
    fprintf(stdout, "Coef result: %f\n", sum);
    
    return result;
}

