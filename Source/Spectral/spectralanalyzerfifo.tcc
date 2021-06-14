#include <cstddef>
#include <math.h>
#include <string.h>

#include "spectral.h"

template <size_t FFT_SIZE, size_t OVERLAP, size_t WINDOW_SIZE>
void SpectralAnalyzerFifo<FFT_SIZE, OVERLAP, WINDOW_SIZE>::Init(
    SPECTRAL_WINDOW window_type,
    size_t          sample_rate)
{
    status_ = STATUS::OK;
    state_  = STATE::INIT;

    float *analwinhalf, *analwinbase;
    float  sum;
    int    halfwinsize;
    int    i, Mf /*,Lf*/;

    unsigned int N       = FFT_SIZE;
    unsigned int M       = WINDOW_SIZE;
    unsigned int overlap = OVERLAP;

    if(overlap < 48 || overlap <= 10) /* 10 is a guess.... */
    {
        // return pvssanalset(csound, p);

        // TODO -- fully implement sliding
        status_ = STATUS::E_SLIDING_NOT_IMPLEMENTED;
        return;

        InitSliding(window_type, sample_rate, 48);
        return;
    }

    if(N & (N - 1))
    {
        status_ = STATUS::E_FFT_NOT_POWER;
        return;
    }

    if(N <= 32)
    {
        status_ = STATUS::E_FFT_TOO_SMALL;
        return;
        // return csound->InitError(csound,
        //                          Str("pvsanal: fftsize of 32 is too small!\n"));
    }

    /* check N for powof2? CARL fft routines and FFTW are not limited to that */
    N = N + N % 2; /* Make N even */

    if(M < N)
    {
        //  csound->Warning(csound,
        //                          Str("pvsanal: window size too small for fftsize"));
        status_ = STATUS::E_WINDOW_TOO_SMALL;
        return;
    }

    if(overlap > N / 2)
    {
        status_ = STATUS::W_OVERLAP_TOO_BIG;
        overlap = (int)N / 2;
        // return csound->InitError(csound,
        //                          Str("pvsanal: overlap too big for fft size\n"));
    }

    half_overlap_  = overlap;
    input_segment_ = overlapbuf_;
    // Initially, the analyzer will be idle, waiting for
    // the input to be filled. When it's actually processing,
    // the processSegment_ and inputSegment_ will be different.
    process_segment_ = input_segment_;

    // NOTE -- only Hamming and Hann are supported
    if(window_type != SPECTRAL_WINDOW::HAMMING
       && window_type != SPECTRAL_WINDOW::HANN)
    {
        status_    = STATUS::W_INVALID_WINDOW;
        window_type = SPECTRAL_WINDOW::HAMMING;
    }

    halfwinsize = M / 2;
    buflen_     = M * 4;
    // arate = (float)(csound->esr / (float) overlap);
    // fund = (float)(csound->esr / (float) N);
    float arate = (float)(sample_rate / (float)overlap);
    // float fund = (float)(sampleRate / (float) N);

    // int nBins = N/2 + 1;

    /* we can exclude/simplify all sorts of stuff in CARL
     * as we will never do time-scaling with this opcode
     */
    /*Lf =*/
    Mf = 1 - M % 2;

    // // NOTE -- take note of these and how they reflect on static allocation
    // csound->AuxAlloc(csound, overlap * sizeof(float), &overlapbuf);
    // csound->AuxAlloc(csound, (N+2) * sizeof(float), &analbuf);
    // csound->AuxAlloc(csound, (M+Mf) * sizeof(float), &analwinbuf);
    // csound->AuxAlloc(csound, nBins * sizeof(float), &oldInPhase);
    // csound->AuxAlloc(csound, buflen_ * sizeof(float), &input);
    // /* the signal itself */
    // csound->AuxAlloc(csound, (N+2) * sizeof(float), &fsig_.frame);

    /* make the analysis window*/
    analwinbase = analwinbuf_;
    analwinhalf = analwinbase + halfwinsize;

    // if (PVS_CreateWindow(csound, analwinhalf, windowType, M) != OK)
    //   return NOTOK;

    SpectralWindow(analwinhalf, window_type, M);

    for(i = 1; i <= halfwinsize; i++)
        *(analwinhalf - i) = *(analwinhalf + i - Mf);
    if(M > N)
    {
        float dN = (float)N;
        /*  sinc function */
        if(Mf)
            *analwinhalf *= (float)(dN * sin(HALFPI_F / dN) / (HALFPI_F));
        for(i = 1; i <= halfwinsize; i++)
            *(analwinhalf + i)
                *= (float)(dN * sin((float)(PI_F * (i + 0.5 * Mf) / dN))
                           / (PI_F * (i + 0.5 * Mf)));
        for(i = 1; i <= halfwinsize; i++)
            *(analwinhalf - i) = *(analwinhalf + i - Mf);
    }
    /* get net amp */
    sum = 0.0f;

    for(i = -halfwinsize; i <= halfwinsize; i++)
        sum += *(analwinhalf + i);
    sum = 2.0f / sum; /* factor of 2 comes in later in trig identity */
    for(i = -halfwinsize; i <= halfwinsize; i++)
        *(analwinhalf + i) *= sum;


    /*    invR = (float)(FL(1.0) / csound->esr); */
    RoverTwoPi_ = (float)(arate / TWOPI_F);
    TwoPioverR_ = (float)(TWOPI_F / arate);
    Fexact_     = (float)(sample_rate / (float)N);
    nI_         = -((int64_t)(halfwinsize / overlap))
          * overlap; /* input time (in samples) */
    /*Dd = halfwinsize + nI_ + 1;                     */
    /* in streaming mode, Dd = ovelap all the time */
    Ii_     = 0;
    IOi_    = 0;
    nextIn_ = input_;
    inptr_  = 0;
    /* and finally, set up the output signal */
    fsig_out_.N          = N;
    fsig_out_.overlap    = overlap;
    fsig_out_.winsize    = M;
    fsig_out_.wintype    = window_type;
    fsig_out_.framecount = 1;
    fsig_out_.format     = SPECTRAL_FORMAT::AMP_FREQ; /* only this, for now */
    fsig_out_.sliding    = false;

    // if (!(N & (N - 1))) /* if pow of two use this */
    //  setup = csound->RealFFT2Setup(csound,N,FFT_FWD);
    // return OK;

    sample_rate_ = sample_rate;
    fft_.Init();
}

template <size_t FFT_SIZE, size_t OVERLAP, size_t WINDOW_SIZE>
void SpectralAnalyzerFifo<FFT_SIZE, OVERLAP, WINDOW_SIZE>::InitSliding(
    SPECTRAL_WINDOW window_type,
    size_t          sample_rate,
    size_t          block)
{
    /* opcode params */
    int N = WINDOW_SIZE;
    int NB;
    int i;

    if(N <= 0)
    {
        status_ = STATUS::E_WINDOW_TOO_SMALL;
        return;
    }
    // if (N<=0) return csound->InitError(csound, Str("Invalid window size"));
    /* deal with iinit and iformat later on! */

    N  = N + N % 2; /* Make N even */
    NB = N / 2 + 1; /* Number of bins */

    /* Need space for NB complex numbers for each of ksmps */

    // NOTE -- static allocation makes sliding particularly difficult
    // if (fsig_.frame==NULL ||
    //     block*(N+2)*sizeof(float) > (uint32_t)fsig_.frame.size)
    //   csound->AuxAlloc(csound, block*(N+2)*sizeof(float),&fsig_.frame);
    // else memset(fsig_.frame, 0, block*(N+2)*sizeof(float));

    /* Space for remembering samples */
    //   if (input.auxp==NULL ||
    //       N*sizeof(float) > (uint32_t)input.size)
    //     csound->AuxAlloc(csound, N*sizeof(float),&input);
    //   else memset(input.auxp, 0, N*sizeof(float));
    //   csound->AuxAlloc(csound, NB * sizeof(float), &oldInPhase);

    //  if (analwinbuf.auxp==NULL ||
    //       NB*sizeof(CMPLX) > (uint32_t)analwinbuf.size)
    //     csound->AuxAlloc(csound, NB*sizeof(CMPLX),&analwinbuf);
    //   else memset(analwinbuf.auxp, 0, NB*sizeof(CMPLX));

    inptr_   = 0; /* Pointer in circular buffer */
    fsig_out_.NB = Ii_ = NB;
    fsig_out_.wintype  = window_type;
    fsig_out_.format   = SPECTRAL_FORMAT::AMP_FREQ; /* only this, for now */
    fsig_out_.N = nI_ = N;
    fsig_out_.sliding = true;

    // NOTE -- avoiding allocation but still filling trig buffers
    /* Need space for NB sines, cosines and a scatch phase area */
    // if (trig.auxp==NULL ||
    //     (2*NB)*sizeof(float) > (uint32_t)trig.size)
    //   csound->AuxAlloc(csound,(2*NB)*sizeof(float),&trig);
    {
        float  dc = cos(TWOPI_F / (float)N);
        float  ds = sin(TWOPI_F / (float)N);
        float *c  = trig_;
        float *s  = c + NB;
        cosine_   = c;
        sine_     = s;
        c[0]      = 1.0;
        s[0]      = 0.0; // assignment to s unnecessary as auxalloc zeros
                         /*
          direct computation of c and s may be better for large n
          c[i] = cos(2*PI*i/n);
          s[i] = sin(2*PI*i/n);
          if (i % 16 == 15) {
          c[i] = cos(2*PI*(i+1)/n);
          s[i] = sin(2*PI*(i+1)/n);
        */
        for(i = 1; i < NB; i++)
        {
            c[i] = dc * c[i - 1] - ds * s[i - 1];
            s[i] = ds * c[i - 1] + dc * s[i - 1];
        }
        /*       for (i=0; i<NB; i++)  */
        /*         printf("c[%d] = %f   \ts[%d] = %f\n", i, c[i], i, s[i]); */
    }
    // return OK;

    sample_rate_ = sample_rate;
    // The fft isn't used in sliding mode
    // fft_.Init();
}

template <size_t FFT_SIZE, size_t OVERLAP, size_t WINDOW_SIZE>
SpectralBuffer<FFT_SIZE>&
SpectralAnalyzerFifo<FFT_SIZE, OVERLAP, WINDOW_SIZE>::Process()
{
    while(state_ != STATE::PROCESSING) {}
    GenerateFrame();
    state_ = STATE::IDLE;
    fsig_out_.framecount++;
    return fsig_out_;


    // // float *ain;
    // // unsigned int offset = h.insdshead->ksmps_offset;
    // // unsigned int early  = h.insdshead->ksmps_no_end;
    // // unsigned int i, nsmps = block;
    // unsigned int early  = 0;
    // unsigned int offset = 0;
    // unsigned int i, nsmps = size;

    // // ain = ain;

    // // NOTE -- keep this in mind if we decide to dynamically allocate
    // // if (input.auxp==NULL) {
    // //   return csound->PerfError(csound,&(h),
    // //                            Str("pvsanal: Not Initialised.\n"));
    // // }
    // // overlap is a pointer to a buffer, so no need to cast?
    // // int over = (int)*overlap;
    // if(fsig_.overlap < (int)nsmps || fsig_.overlap < 10)
    // {
    //     ProcessSliding(in, size);
    // }
    // else
    // {
    //     nsmps -= early;

    //     for(i = offset; i < nsmps; i++)
    //         Tick(in[i]);
    // }

    // return fsig_;
}

template <size_t FFT_SIZE, size_t OVERLAP, size_t WINDOW_SIZE>
void SpectralAnalyzerFifo<FFT_SIZE, OVERLAP, WINDOW_SIZE>::Sample(float sample)
{
    if(inptr_ == fsig_out_.overlap)
    {
        inptr_ = 0;
        switch(state_)
        {
            case STATE::INIT:
                input_segment_ = overlapbuf_ + half_overlap_;
                state_        = STATE::PROCESSING;
                break;
            case STATE::IDLE:
                input_segment_ = (input_segment_ == overlapbuf_)
                                    ? overlapbuf_ + half_overlap_
                                    : overlapbuf_;
                process_segment_ = (process_segment_ == overlapbuf_)
                                      ? overlapbuf_ + half_overlap_
                                      : overlapbuf_;
                state_ = STATE::PROCESSING;
                break;
            case STATE::PROCESSING: status_ = STATUS::W_BUFFER_UNDERFLOW; break;
            default: status_ = STATUS::W_INVALID_STATE; break;
        }
    }

    input_segment_[inptr_++] = sample;
}

// template <size_t FFT_SIZE, size_t OVERLAP, size_t WINDOW_SIZE>
// void SpectralAnalyzerFifo<FFT_SIZE, OVERLAP, WINDOW_SIZE>::Tick(float sample)
// {
//     if(inptr_ == fsig_out_.overlap)
//     {
//         GenerateFrame();
//         fsig_out_.framecount++;
//         inptr_ = 0;
//     }
//     //printf("inptr_ = %d fsig_.overlap=%d\n", inptr_, fsig_.overlap);
//     overlapbuf_[inptr_++] = sample;
// }

template <size_t FFT_SIZE, size_t OVERLAP, size_t WINDOW_SIZE>
void SpectralAnalyzerFifo<FFT_SIZE, OVERLAP, WINDOW_SIZE>::ProcessSliding(
    const float *in,
    size_t       size)
{
    // float *ain;
    int      NB   = Ii_, loc;
    int      N    = fsig_out_.N;
    float *  data = input_;
    Complex *fw   = (Complex *)
        analwinbuf_; // casting this regular float buffer to complex values is annoying
    float *c = cosine_;
    float *s = sine_;
    float *h = oldInPhase_;

    unsigned int offset = 0;
    unsigned int early  = 0;
    unsigned int i, nsmps = size;

    // NOTE -- keep this in mind if we decide to dynamically allocate
    // if (data==NULL) {
    //   return csound->PerfError(csound,&(h),
    //                            Str("pvsanal: Not Initialised.\n"));
    // }

    // ain = ain;               /* The input samples */
    loc = inptr_; /* Circular buffer */
    nsmps -= early;
    for(i = offset; i < nsmps; i++)
    {
        float    re, im, dx;
        Complex *ff;
        int      j;

        /*       printf("%d: in = %f\n", i, *ain); */
        // dx = *ain - data[loc];    /* Change in sample */
        // data[loc] = *ain++;       /* Remember input sample */
        dx        = *in - data[loc];
        data[loc] = *in++;
        /* get the frame for this sample */

        ff = (Complex *)(fsig_out_.frame) + i * NB;
        /* fw is the current frame at this sample */
        for(j = 0; j < NB; j++)
        {
            float ci = c[j], si = s[j];
            re      = fw[j].a + dx;
            im      = fw[j].b;
            fw[j].a = ci * re - si * im;
            fw[j].b = ci * im + si * re;
        }
        loc++;
        if(loc == nI_)
            loc = 0; /* Circular buffer */
        /* apply window and transfer to ff buffer*/
        /* Rectang :Fw_t =     F_t                          */
        /* Hamming :Fw_t = 0.54F_t - 0.23[ F_{t-1}+F_{t+1}] */
        /* Hamming :Fw_t = 0.5 F_t - 0.25[ F_{t-1}+F_{t+1}] */
        /* Blackman:Fw_t = 0.42F_t - 0.25[ F_{t-1}+F_{t+1}]+0.04[F_{t-2}+F_{t+2}] */
        /* Blackman_exact:Fw_t = 0.42659071367153912296F_t
         - 0.24828030954428202923 [F_{t-1}+F_{t+1}]
         + 0.038424333619948409286 [F_{t-2}+F_{t+2}]      */
        /* Nuttall_C3:Fw_t = 0.375  F_t - 0.25[ F_{t-1}+F_{t+1}] +
                                      0.0625 [F_{t-2}+F_{t+2}] */
        /* BHarris_3:Fw_t = 0.44959 F_t - 0.24682[ F_{t-1}+F_{t+1}] +
                                      0.02838 [F_{t-2}+F_{t+2}] */
        /* BHarris_min:Fw_t = 0.42323 F_t - 0.2486703 [ F_{t-1}+F_{t+1}] +
                                      0.0391396 [F_{t-2}+F_{t+2}] */
        switch(fsig_out_.wintype)
        {
            case SPECTRAL_WINDOW::HAMMING:
                for(j = 0; j < NB; j++)
                {
                    ff[j].a = 0.54f * fw[j].a;
                    ff[j].b = 0.54f * fw[j].b;
                }
                for(j = 1; j < NB - 1; j++)
                {
                    ff[j].a -= 0.23f * (fw[j + 1].a + fw[j - 1].a);
                    ff[j].b -= 0.23f * (fw[j + 1].b + fw[j - 1].b);
                }
                ff[0].a -= 0.46f * fw[1].a;
                ff[NB - 1].a -= 0.46f * fw[NB - 2].a;
                break;
            case SPECTRAL_WINDOW::HANN:
                for(j = 0; j < NB; j++)
                {
                    ff[j].a = 0.5f * fw[j].a;
                    ff[j].b = 0.5f * fw[j].b;
                }
                for(j = 1; j < NB - 1; j++)
                {
                    ff[j].a -= 0.25f * (fw[j + 1].a + fw[j - 1].a);
                    ff[j].b -= 0.25f * (fw[j + 1].b + fw[j - 1].b);
                }
                ff[0].a -= 0.5f * fw[1].a;
                ff[NB - 1].a -= 0.5f * fw[NB - 2].a;
                break;
            default:
                status_ = STATUS::W_INVALID_WINDOW;
                // csound->Warning(csound,
                //                 Str("Unknown window type; replaced by rectangular\n"));
                /* FALLTHRU */
            case SPECTRAL_WINDOW::RECT:
                memcpy(ff, fw, NB * sizeof(Complex));
                /* for (j=0; j<NB; j++) { */
                /*   ff[j].a = fw[j].a; */
                /*   ff[j].b = fw[j].b; */
                /* } */
                break;
            case SPECTRAL_WINDOW::BLACKMAN:
                for(j = 0; j < NB; j++)
                {
                    ff[j].a = 0.42f * fw[j].a;
                    ff[j].b = 0.42f * fw[j].b;
                }
                for(j = 1; j < NB - 1; j++)
                {
                    ff[j].a -= 0.25f * (fw[j + 1].a + fw[j - 1].a);
                    ff[j].b -= 0.25f * (fw[j + 1].b + fw[j - 1].b);
                }
                for(j = 2; j < NB - 2; j++)
                {
                    ff[j].a += 0.04f * (fw[j + 2].a + fw[j - 2].a);
                    ff[j].b += 0.04f * (fw[j + 2].b + fw[j - 2].b);
                }
                ff[0].a += -0.5f * fw[1].a + 0.08f * fw[2].a;
                ff[NB - 1].a += -0.5f * fw[NB - 2].a + 0.08f * fw[NB - 3].a;
                ff[1].a += -0.5f * fw[2].a + 0.08f * fw[3].a;
                ff[NB - 2].a += -0.5f * fw[NB - 3].a + 0.08f * fw[NB - 4].a;
                break;
            case SPECTRAL_WINDOW::BLACKMAN_EXACT:
                for(j = 0; j < NB; j++)
                {
                    ff[j].a = 0.42659071367153912296f * fw[j].a;
                    ff[j].b = 0.42659071367153912296f * fw[j].b;
                }
                for(j = 1; j < NB - 1; j++)
                {
                    ff[j].a -= 0.49656061908856405847f * 0.5f
                               * (fw[j + 1].a + fw[j - 1].a);
                    ff[j].b -= 0.49656061908856405847f * 0.5f
                               * (fw[j + 1].b + fw[j - 1].b);
                }
                for(j = 2; j < NB - 2; j++)
                {
                    ff[j].a += 0.076848667239896818573f * 0.5f
                               * (fw[j + 2].a + fw[j - 2].a);
                    ff[j].b += 0.076848667239896818573f * 0.5f
                               * (fw[j + 2].b + fw[j - 2].b);
                }
                ff[0].a += -0.49656061908856405847f * fw[1].a
                           + 0.076848667239896818573f * fw[2].a;
                ff[NB - 1].a += -0.49656061908856405847f * fw[NB - 2].a
                                + 0.076848667239896818573f * fw[NB - 3].a;
                ff[1].a += -0.49656061908856405847f * fw[2].a
                           + 0.076848667239896818573f * fw[3].a;
                ff[NB - 2].a += -0.49656061908856405847f * fw[NB - 3].a
                                + 0.076848667239896818573f * fw[NB - 4].a;
                break;
            case SPECTRAL_WINDOW::NUTTALLC3:
                for(j = 0; j < NB; j++)
                {
                    ff[j].a = 0.375f * fw[j].a;
                    ff[j].b = 0.375f * fw[j].b;
                }
                for(j = 1; j < NB - 1; j++)
                {
                    ff[j].a -= 0.5f * 0.5f * (fw[j + 1].a + fw[j - 1].a);
                    ff[j].b -= 0.5f * 0.5f * (fw[j + 1].b + fw[j - 1].b);
                }
                for(j = 2; j < NB - 2; j++)
                {
                    ff[j].a += 0.125f * 0.5f * (fw[j + 2].a + fw[j - 2].a);
                    ff[j].b += 0.125f * 0.5f * (fw[j + 2].b + fw[j - 2].b);
                }
                ff[0].a += -0.5f * fw[1].a + 0.125f * fw[2].a;
                ff[NB - 1].a += -0.5f * fw[NB - 2].a + 0.125f * fw[NB - 3].a;
                ff[1].a += -0.5f * fw[2].a + 0.125f * fw[3].a;
                ff[NB - 2].a += -0.5f * fw[NB - 3].a + 0.125f * fw[NB - 4].a;
                ff[1].a = 0.5 * (fw[2].a + fw[0].a); /* HACK???? */
                ff[1].b = 0.5 * (fw[2].b + fw[0].b);
                break;
            case SPECTRAL_WINDOW::BHARRIS_3:
                for(j = 0; j < NB; j++)
                {
                    ff[j].a = 0.44959f * fw[j].a;
                    ff[j].b = 0.44959f * fw[j].b;
                }
                for(j = 1; j < NB - 1; j++)
                {
                    ff[j].a -= 0.49364f * 0.5f * (fw[j + 1].a + fw[j - 1].a);
                    ff[j].b -= 0.49364f * 0.5f * (fw[j + 1].b + fw[j - 1].b);
                }
                for(j = 2; j < NB - 2; j++)
                {
                    ff[j].a += 0.05677f * 0.5f * (fw[j + 2].a + fw[j - 2].a);
                    ff[j].b += 0.05677f * 0.5f * (fw[j + 2].b + fw[j - 2].b);
                }
                ff[0].a += -0.49364f * fw[1].a + 0.05677f * fw[2].a;
                ff[NB - 1].a
                    += -0.49364f * fw[NB - 2].a + 0.05677f * fw[NB - 3].a;
                ff[1].a += -0.49364f * fw[2].a + 0.05677f * fw[3].a;
                ff[NB - 2].a
                    += -0.49364f * fw[NB - 3].a + 0.05677f * fw[NB - 4].a;
                ff[1].a = 0.5 * (fw[2].a + fw[0].a); /* HACK???? */
                ff[1].b = 0.5 * (fw[2].b + fw[0].b);
                break;
            case SPECTRAL_WINDOW::BHARRIS_MIN:
                for(j = 0; j < NB; j++)
                {
                    ff[j].a = 0.42323f * fw[j].a;
                    ff[j].b = 0.42323f * fw[j].b;
                }
                for(j = 1; j < NB - 1; j++)
                {
                    ff[j].a -= 0.4973406f * 0.5f * (fw[j + 1].a + fw[j - 1].a);
                    ff[j].b -= 0.4973406f * 0.5f * (fw[j + 1].b + fw[j - 1].b);
                }
                for(j = 2; j < NB - 2; j++)
                {
                    ff[j].a += 0.0782793f * 0.5f * (fw[j + 2].a + fw[j - 2].a);
                    ff[j].b += 0.0782793f * 0.5f * (fw[j + 2].b + fw[j - 2].b);
                }
                ff[0].a += -0.4973406f * fw[1].a + 0.0782793f * fw[2].a;
                ff[NB - 1].a
                    += -0.4973406f * fw[NB - 2].a + 0.0782793f * fw[NB - 3].a;
                ff[1].a += -0.4973406f * fw[2].a + 0.0782793f * fw[3].a;
                ff[NB - 2].a
                    += -0.4973406f * fw[NB - 3].a + 0.0782793f * fw[NB - 4].a;
                ff[1].a = 0.5 * (fw[2].a + fw[0].a); /* HACK???? */
                ff[1].b = 0.5 * (fw[2].b + fw[0].b);
                break;
        }
        /*       if (i==9) { */
        /*         printf("Frame as Amp/Freq %d\n", i); */
        /*         for (j = 0; j < NB; j++) */
        /*           printf("%d: %f\t%f\n", j, ff[j].a, ff[j].b); */
        /*       } */
        for(j = 0; j < NB; j++)
        { /* Convert to AMP_FREQ */
            float thismag  = hypot(ff[j].a, ff[j].b);
            float phase    = atan2(ff[j].b, ff[j].a);
            float angleDif = phase - h[j];
            h[j]           = phase;
            /*subtract expected phase difference */
            angleDif -= (float)j * TWOPI_F / N;
            angleDif = mod2Pi(angleDif);
            angleDif = angleDif * N / TWOPI_F;
            ff[j].a  = thismag;

            // ff[j].b = csound->esr * (j + angleDif)/N;
            ff[j].b = sample_rate_ * (j + angleDif) / N;
        }
        /*       if (i==9) { */
        /*         printf("Frame as Amp/Freq %d\n", i); */
        /*         for (j = 0; j < NB; j++) */
        /*           printf("%d: %f\t%f\n", j, ff[j].a, ff[j].b); */
        /*       } */
    }

    inptr_ = loc;
    // return OK;
}

template <size_t FFT_SIZE, size_t OVERLAP, size_t WINDOW_SIZE>
void SpectralAnalyzerFifo<FFT_SIZE, OVERLAP, WINDOW_SIZE>::GenerateFrame()
{
    int    got, tocp, i, j, k, ii;
    int    N          = fsig_out_.N;
    int    N2         = N / 2;
    int    analWinLen = fsig_out_.winsize / 2;
    int    synWinLen  = analWinLen;
    float *ofp; /* RWD MUST be 32bit */
    float *fp;
    float *anal           = analbuf_;
    float *analOut        = analbufOut_;
    float *tempInput      = input_;
    float *analWindow     = analwinbuf_ + analWinLen;
    float *tempOldInPhase = oldInPhase_;
    float  angleDif, real, imag, phase;
    float  rratio;

    got = fsig_out_.overlap; /*always assume */
    // fp   = overlapbuf_;
    fp   = process_segment_;
    tocp = (got <= tempInput + buflen_ - nextIn_
                ? got
                : tempInput + buflen_ - nextIn_);
    got -= tocp;
    while(tocp-- > 0)
        *(nextIn_++) = *fp++;

    if(got > 0)
    {
        nextIn_ -= buflen_;
        while(got-- > 0)
            *nextIn_++ = *fp++;
    }
    if(nextIn_ >= (tempInput + buflen_))
        nextIn_ -= buflen_;

    /* analysis: The analysis subroutine computes the complex output at
       time n of (N/2 + 1) of the phase vocoder channels.  It operates
       on input samples (n - analWinLen) thru (n + analWinLen) and
       expects to find these in input[(n +- analWinLen) mod ibuflen].
       It expects analWindow to point to the center of a
       symmetric window of length (2 * analWinLen +1).  It is the
       responsibility of the main program to ensure that these values
       are correct!  The results are returned in anal as succesive
       pairs of real and imaginary values for the lowest (N/2 + 1)
       channels.   The subroutines fft and reals together implement
       one efficient FFT call for a real input sequence.  */

    /* for (i = 0; i < N+2; i++)
     *(anal + i) = FL(0.0);  */     /*initialize*/
    memset(anal, 0, sizeof(float) * (N + 2));

    j = (nI_ - analWinLen - 1 + buflen_) % buflen_; /*input pntr*/

    k = nI_ - analWinLen - 1; /*time shift*/
    while(k < 0)
        k += N;
    k = k % N;
    for(i = -analWinLen; i <= analWinLen; i++)
    {
        if(++j >= buflen_)
            j -= buflen_;
        if(++k >= N)
            k -= N;
        /* *(anal + k) += *(analWindow + i) * *(input + j); */
        anal[k] += analWindow[i] * tempInput[j];
    }
    // if (!(N & (N - 1))) { // NOTE -- if N is a power of two (including 2^0)
    //   /* csound->RealFFT(csound, anal, N);*/
    //   csound->RealFFT2(csound,setup,anal);
    //   anal[N] = anal[1];
    //   anal[1] = anal[N + 1] = 0.0f;
    // }
    // else
    //   csound->RealFFTnp2(csound, anal, N);

    //////////////////////////////////////////////////////////
    // Custom FFT section
    //////////////////////////////////////////////////////////

    if(N != FFT_SIZE)
    {
        int num_passes = GetPasses(N);
        fft_.Direct(anal, analOut, num_passes);
    }
    else
    {
        fft_.Direct(anal, analOut);
    }
    Interlace(analOut, anal, N);

    //////////////////////////////////////////////////////////
    // Custom FFT section end
    //////////////////////////////////////////////////////////

    /* conversion: The real and imaginary values in anal are converted to
       magnitude and angle-difference-per-second (assuming an
       intermediate sampling rate of rIn) and are returned in
       anal. */

    /*if (format==PVS_AMP_FREQ) {*/
    for(i = ii = 0 /*,i0=anal,i1=anal+1,oi=tempOldInPhase*/; i <= N2;
        i++, ii += 2 /*i0+=2,i1+=2, oi++*/)
    {
        real             = anal[ii] /* *i0 */;
        imag             = anal[ii + 1] /* *i1 */;
        /**i0*/ anal[ii] = hypotf(real, imag);
        /* phase unwrapping */
        /*if (*i0 == 0.)*/
        if(/* *i0 */ anal[ii] < 1.0E-10f)
            angleDif = 0.0f;
        else
        {
            rratio   = atan2((float)imag, (float)real);
            angleDif = (phase = (float)rratio) - /**oi*/ tempOldInPhase[i];
            /* *oi */ tempOldInPhase[i] = phase;
        }

        if(angleDif > PI_F)
            angleDif = angleDif - TWOPI_F;
        if(angleDif < -PI_F)
            angleDif = angleDif + TWOPI_F;

        /* add in filter center freq.*/
        /* *i1 */ anal[ii + 1] = angleDif * RoverTwoPi_ + ((float)i * Fexact_);
    }
    /* } */
    /* else must be PVOC_COMPLEX */
    fp  = anal;
    ofp = fsig_out_.frame; /* RWD MUST be 32bit */
    for(i = 0; i < N + 2; i++)
        /* *ofp++ = (float)(*fp++); */
        ofp[i] = (float)fp[i];

    nI_ += fsig_out_.overlap; /* increment time */
    if(nI_ > (synWinLen + fsig_out_.overlap))
        Ii_ = /*I*/ fsig_out_.overlap;
    else if(nI_ > synWinLen)
        Ii_ = nI_ - synWinLen;
    else
    {
        Ii_ = 0;

        /*  for (i=nO+synWinLen; i<buflen_; i++)
            if (i > 0)
            *(output+i) = 0.0f;
            */
    }

    IOi_ = Ii_;
}
