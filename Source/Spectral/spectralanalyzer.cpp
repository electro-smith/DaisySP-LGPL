#include <math.h>
#include <string.h>
#include "dsp.h"
#include "spectral.h"
#include "spectralanalyzer.h"
// #include "shy_fft.h"

using namespace daicsp;

void SpectralAnalyzer::Init(
  int fftsize, 
  int overlap, 
  int windowSize, 
  SPECTRAL_WINDOW windowType, 
  size_t sampleRate,
  size_t blockSize
)
{

    float *analwinhalf,*analwinbase;
    float sum;
    int halfwinsize;
    int i,Mf/*,Lf*/;

    unsigned int N = fftsize;
    unsigned int M = windowSize;

    if (overlap < (int) blockSize || overlap <= 10) /* 10 is a guess.... */
    {
      // return pvssanalset(csound, p);
      // TODO -- put in extra special init here
      InitSmall();
      return;
    } 

    // TODO -- error handling
    if (N <= 32)
    {
      // return csound->InitError(csound,
      //                          Str("pvsanal: fftsize of 32 is too small!\n"));
    }
      

    /* check N for powof2? CARL fft routines and FFTW are not limited to that */
    N = N  + N%2;       /* Make N even */
    // TODO -- error handling
    if (M < N) 
    {
      //  csound->Warning(csound,
      //                          Str("pvsanal: window size too small for fftsize"));
       M = N;
    }

    // TODO -- error handling
    if (overlap > (int) N / 2)
    {
      // return csound->InitError(csound,
      //                          Str("pvsanal: overlap too big for fft size\n"));
    }

    halfwinsize = M/2;
    buflen = M*4;
    // arate = (float)(csound->esr / (float) overlap);
    // fund = (float)(csound->esr / (float) N);
    float arate = (float)(sampleRate / (float) overlap);
    // float fund = (float)(sampleRate / (float) N);

    // int nBins = N/2 + 1;

    /* we can exclude/simplify all sorts of stuff in CARL
     * as we will never do time-scaling with this opcode
     */
    /*Lf =*/ 
    Mf = 1 - M%2;

    // // NOTE -- take note of these and how they reflect on static allocation
    // csound->AuxAlloc(csound, overlap * sizeof(float), &overlapbuf);
    // csound->AuxAlloc(csound, (N+2) * sizeof(float), &analbuf);
    // csound->AuxAlloc(csound, (M+Mf) * sizeof(float), &analwinbuf);
    // csound->AuxAlloc(csound, nBins * sizeof(float), &oldInPhase);
    // csound->AuxAlloc(csound, buflen * sizeof(float), &input);
    // /* the signal itself */
    // csound->AuxAlloc(csound, (N+2) * sizeof(float), &fsig_.frame);

    /* make the analysis window*/
    analwinbase = analwinbuf;
    analwinhalf = analwinbase + halfwinsize;

    // TODO -- error handling
    // if (PVS_CreateWindow(csound, analwinhalf, windowType, M) != OK)
    //   return NOTOK;

    // TODO -- manage the window creation
    SpectralWindow(analwinhalf, windowType, M);

    for (i = 1; i <= halfwinsize; i++)
      *(analwinhalf - i) = *(analwinhalf + i - Mf);
    if (M > N) {
      double dN = (double)N;
      /*  sinc function */
      if (Mf)
        *analwinhalf *= (float)(dN * sin(HALFPI_F/dN) / (HALFPI_F));
      for (i = 1; i <= halfwinsize; i++)
        *(analwinhalf + i) *= (float)
          (dN * sin((double)(PI_F*(i+0.5*Mf)/dN)) / (PI_F*(i+0.5*Mf)));
      for (i = 1; i <= halfwinsize; i++)
        *(analwinhalf - i) = *(analwinhalf + i - Mf);
    }
    /* get net amp */
    sum = 0.0f;

    for (i = -halfwinsize; i <= halfwinsize; i++)
      sum += *(analwinhalf + i);
    sum = 2.0f / sum;  /* factor of 2 comes in later in trig identity */
    for (i = -halfwinsize; i <= halfwinsize; i++)
      *(analwinhalf + i) *= sum;


  /*    invR = (float)(FL(1.0) / csound->esr); */
    RoverTwoPi = (float)(arate / TWOPI_F);
    TwoPioverR = (float)(TWOPI_F / arate);
    Fexact =  (float)(sampleRate / (float)N);
    nI = -((int64_t)(halfwinsize/overlap))*overlap; /* input time (in samples) */
    /*Dd = halfwinsize + nI + 1;                     */
    /* in streaming mode, Dd = ovelap all the time */
    Ii = 0;
    IOi = 0;
    nextIn = input;
    inptr = 0;
    /* and finally, set up the output signal */
    fsig_.N = N;
    fsig_.overlap = overlap;
    fsig_.winsize = M;
    fsig_.wintype = windowType;
    fsig_.framecount = 1;
    fsig_.format = SPECTRAL_FORMAT::AMP_FREQ;      /* only this, for now */
    fsig_.sliding = 0;

    // if (!(N & (N - 1))) /* if pow of two use this */
    //  setup = csound->RealFFT2Setup(csound,N,FFT_FWD);
    // return OK;

    sr_ = sampleRate;

    fft_.Init();
}

void SpectralAnalyzer::InitSmall() 
{
  
}

SpectralBuffer& SpectralAnalyzer::Process(const float* in, size_t size)
{

    // float *ain;
    // unsigned int offset = h.insdshead->ksmps_offset;
    // unsigned int early  = h.insdshead->ksmps_no_end;
    // unsigned int i, nsmps = CS_KSMPS;
    unsigned int early = 0;
    unsigned int offset = 0;
    unsigned int i, nsmps = size;

    // ain = ain;

    // TODO -- error reporting
    // if (input.auxp==NULL) {
    //   return csound->PerfError(csound,&(h),
    //                            Str("pvsanal: Not Initialised.\n"));
    // }
    // overlap is a pointer to a buffer, so no need to cast?
    // int over = (int)*overlap;
    // This seems backwards -- shouldn't it tick when it's still collecting samples?
    if (fsig_.overlap<(int)nsmps || fsig_.overlap<10) {
      Analyze(in, size);
    } else {
      nsmps -= early;

      for (i=offset; i < nsmps; i++)
        // Tick(ain[i]);
        Tick(in[i]);
    }

    return fsig_;
}

void SpectralAnalyzer::Analyze(const float* in, size_t size)
{
    // float *ain;
    int NB = Ii, loc;
    int N = fsig_.N;
    float *data = input;
    Complex *fw = (Complex*) analwinbuf; // casting this regular float buffer to complex values is annoying
    float *c = cosine;
    float *s = sine;
    float *h = oldInPhase;

    // TODO -- make these csound values work with daisy syntax
    // unsigned int offset = h.insdshead->ksmps_offset;
    // unsigned int early  = h.insdshead->ksmps_no_end;

    unsigned int offset = 0;
    unsigned int early = 0;
    // unsigned int i, nsmps = CS_KSMPS;
    unsigned int i, nsmps = size;

    // TODO -- shall we include error handling?
    // if (data==NULL) {
    //   return csound->PerfError(csound,&(h),
    //                            Str("pvsanal: Not Initialised.\n"));
    // }

    // ain = ain;               /* The input samples */
    loc = inptr;             /* Circular buffer */
    nsmps -= early;
    for (i=offset; i < nsmps; i++) {
      float re, im, dx;
      Complex* ff;
      int j;

/*       printf("%d: in = %f\n", i, *ain); */
      // dx = *ain - data[loc];    /* Change in sample */
      // data[loc] = *ain++;       /* Remember input sample */
      dx = *in - data[loc];
      data[loc] = *in++;
      /* get the frame for this sample */

      // TODO -- is this always a set of complex numbers? if so, just declare it as complex, no?
      ff = (Complex*)(fsig_.frame) + i*NB;
      /* fw is the current frame at this sample */
      for (j = 0; j < NB; j++) {
        float ci = c[j], si = s[j];
        re = fw[j].a + dx;
        im = fw[j].b;
        fw[j].a = ci*re - si*im;
        fw[j].b = ci*im + si*re;
      }
      loc++; if (loc==nI) loc = 0; /* Circular buffer */
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
      switch (fsig_.wintype) {
      case SPECTRAL_WINDOW::HAMMING:
        for (j=0; j<NB; j++) {
          ff[j].a = 0.54f*fw[j].a;
          ff[j].b = 0.54f*fw[j].b;
        }
        for (j=1; j<NB-1; j++) {
          ff[j].a -= 0.23f*(fw[j+1].a + fw[j-1].a);
          ff[j].b -= 0.23f*(fw[j+1].b + fw[j-1].b);
        }
        ff[0].a -= 0.46f*fw[1].a;
        ff[NB-1].a -= 0.46f*fw[NB-2].a;
        break;
      case SPECTRAL_WINDOW::HANN:
        for (j=0; j<NB; j++) {
          ff[j].a = 0.5f*fw[j].a;
          ff[j].b = 0.5f*fw[j].b;
        }
        for (j=1; j<NB-1; j++) {
          ff[j].a -= 0.25f*(fw[j+1].a + fw[j-1].a);
          ff[j].b -= 0.25f*(fw[j+1].b + fw[j-1].b);
        }
        ff[0].a -= 0.5f*fw[1].a;
        ff[NB-1].a -= 0.5f*fw[NB-2].a;
        break;
      default:
        // TODO -- error handling?
        // csound->Warning(csound,
        //                 Str("Unknown window type; replaced by rectangular\n"));
        /* FALLTHRU */
      case SPECTRAL_WINDOW::RECT:
        memcpy(ff, fw, NB*sizeof(Complex));
        /* for (j=0; j<NB; j++) { */
        /*   ff[j].a = fw[j].a; */
        /*   ff[j].b = fw[j].b; */
        /* } */
        break;
      case SPECTRAL_WINDOW::BLACKMAN:
        for (j=0; j<NB; j++) {
          ff[j].a = 0.42f*fw[j].a;
          ff[j].b = 0.42f*fw[j].b;
        }
        for (j=1; j<NB-1; j++) {
          ff[j].a -= 0.25f*(fw[j+1].a + fw[j-1].a);
          ff[j].b -= 0.25f*(fw[j+1].b + fw[j-1].b);
        }
        for (j=2; j<NB-2; j++) {
          ff[j].a += 0.04f*(fw[j+2].a + fw[j-2].a);
          ff[j].b += 0.04f*(fw[j+2].b + fw[j-2].b);
        }
        ff[0].a    += -0.5f*fw[1].a + 0.08f*fw[2].a;
        ff[NB-1].a += -0.5f*fw[NB-2].a + 0.08f*fw[NB-3].a;
        ff[1].a    += -0.5f*fw[2].a + 0.08f*fw[3].a;
        ff[NB-2].a += -0.5f*fw[NB-3].a + 0.08f*fw[NB-4].a;
        break;
    case SPECTRAL_WINDOW::BLACKMAN_EXACT:
        for (j=0; j<NB; j++) {
          ff[j].a = 0.42659071367153912296f*fw[j].a;
          ff[j].b = 0.42659071367153912296f*fw[j].b;
        }
        for (j=1; j<NB-1; j++) {
          ff[j].a -= 0.49656061908856405847f*0.5f*(fw[j+1].a + fw[j-1].a);
          ff[j].b -= 0.49656061908856405847f*0.5f*(fw[j+1].b + fw[j-1].b);
        }
        for (j=2; j<NB-2; j++) {
          ff[j].a += 0.076848667239896818573f*0.5f*(fw[j+2].a + fw[j-2].a);
          ff[j].b += 0.076848667239896818573f*0.5f*(fw[j+2].b + fw[j-2].b);
        }
        ff[0].a    += -0.49656061908856405847f * fw[1].a
                      + 0.076848667239896818573f * fw[2].a;
        ff[NB-1].a += -0.49656061908856405847f * fw[NB-2].a
                      + 0.076848667239896818573f * fw[NB-3].a;
        ff[1].a    += -0.49656061908856405847f * fw[2].a
                      + 0.076848667239896818573f * fw[3].a;
        ff[NB-2].a += -0.49656061908856405847f * fw[NB-3].a
                      + 0.076848667239896818573f * fw[NB-4].a;
        break;
      case SPECTRAL_WINDOW::NUTTALLC3:
        for (j=0; j<NB; j++) {
          ff[j].a = 0.375f*fw[j].a;
          ff[j].b = 0.375f*fw[j].b;
        }
        for (j=1; j<NB-1; j++) {
          ff[j].a -= 0.5f*0.5f*(fw[j+1].a + fw[j-1].a);
          ff[j].b -= 0.5f*0.5f*(fw[j+1].b + fw[j-1].b);
        }
        for (j=2; j<NB-2; j++) {
          ff[j].a += 0.125f*0.5f*(fw[j+2].a + fw[j-2].a);
          ff[j].b += 0.125f*0.5f*(fw[j+2].b + fw[j-2].b);
        }
        ff[0].a    += -0.5f * fw[1].a    + 0.125f * fw[2].a;
        ff[NB-1].a += -0.5f * fw[NB-2].a + 0.125f * fw[NB-3].a;
        ff[1].a    += -0.5f * fw[2].a    + 0.125f * fw[3].a;
        ff[NB-2].a += -0.5f * fw[NB-3].a + 0.125f * fw[NB-4].a;
        ff[1].a = 0.5 * (fw[2].a + fw[0].a); /* HACK???? */
        ff[1].b = 0.5 * (fw[2].b + fw[0].b);
        break;
      case SPECTRAL_WINDOW::BHARRIS_3:
        for (j=0; j<NB; j++) {
          ff[j].a = 0.44959f*fw[j].a;
          ff[j].b = 0.44959f*fw[j].b;
        }
        for (j=1; j<NB-1; j++) {
          ff[j].a -= 0.49364f*0.5f*(fw[j+1].a + fw[j-1].a);
          ff[j].b -= 0.49364f*0.5f*(fw[j+1].b + fw[j-1].b);
        }
        for (j=2; j<NB-2; j++) {
          ff[j].a += 0.05677f*0.5f*(fw[j+2].a + fw[j-2].a);
          ff[j].b += 0.05677f*0.5f*(fw[j+2].b + fw[j-2].b);
        }
        ff[0].a    += -0.49364f * fw[1].a    + 0.05677f * fw[2].a;
        ff[NB-1].a += -0.49364f * fw[NB-2].a + 0.05677f * fw[NB-3].a;
        ff[1].a    += -0.49364f * fw[2].a    + 0.05677f * fw[3].a;
        ff[NB-2].a += -0.49364f * fw[NB-3].a + 0.05677f * fw[NB-4].a;
        ff[1].a = 0.5 * (fw[2].a + fw[0].a); /* HACK???? */
        ff[1].b = 0.5 * (fw[2].b + fw[0].b);
        break;
      case SPECTRAL_WINDOW::BHARRIS_MIN:
        for (j=0; j<NB; j++) {
          ff[j].a = 0.42323f*fw[j].a;
          ff[j].b = 0.42323f*fw[j].b;
        }
        for (j=1; j<NB-1; j++) {
          ff[j].a -= 0.4973406f*0.5f*(fw[j+1].a + fw[j-1].a);
          ff[j].b -= 0.4973406f*0.5f*(fw[j+1].b + fw[j-1].b);
        }
        for (j=2; j<NB-2; j++) {
          ff[j].a += 0.0782793f*0.5f*(fw[j+2].a + fw[j-2].a);
          ff[j].b += 0.0782793f*0.5f*(fw[j+2].b + fw[j-2].b);
        }
        ff[0].a    += -0.4973406f * fw[1].a    + 0.0782793f * fw[2].a;
        ff[NB-1].a += -0.4973406f * fw[NB-2].a + 0.0782793f * fw[NB-3].a;
        ff[1].a    += -0.4973406f * fw[2].a    + 0.0782793f * fw[3].a;
        ff[NB-2].a += -0.4973406f * fw[NB-3].a + 0.0782793f * fw[NB-4].a;
        ff[1].a = 0.5 * (fw[2].a + fw[0].a); /* HACK???? */
        ff[1].b = 0.5 * (fw[2].b + fw[0].b);
        break;
      }
/*       if (i==9) { */
/*         printf("Frame as Amp/Freq %d\n", i); */
/*         for (j = 0; j < NB; j++) */
/*           printf("%d: %f\t%f\n", j, ff[j].a, ff[j].b); */
/*       } */
      for (j = 0; j < NB; j++) { /* Convert to AMP_FREQ */
        float thismag = hypot(ff[j].a, ff[j].b);
        float phase = atan2(ff[j].b, ff[j].a);
        float angleDif  = phase -  h[j];
        h[j] = phase;
            /*subtract expected phase difference */
        angleDif -= (float)j * TWOPI_F/N;
        angleDif =  mod2Pi(angleDif);
        angleDif =  angleDif * N /TWOPI_F;
        ff[j].a = thismag;
        
        // TODO -- not sure what this esr property represents
        // I'll just assume it's sample rate for now
        // ff[j].b = csound->esr * (j + angleDif)/N;
        ff[j].b = sr_ * (j + angleDif)/N;
      }
/*       if (i==9) { */
/*         printf("Frame as Amp/Freq %d\n", i); */
/*         for (j = 0; j < NB; j++) */
/*           printf("%d: %f\t%f\n", j, ff[j].a, ff[j].b); */
/*       } */
    }

    inptr = loc;
    // return OK;
}

void SpectralAnalyzer::Tick(float sample) 
{
    if (inptr == fsig_.overlap) {
      UpdateFrame();
      fsig_.framecount++;
      inptr = 0;
    }
    //printf("inptr = %d fsig_.overlap=%d\n", inptr, fsig_.overlap);
    overlapbuf[inptr++] = sample;
}

void SpectralAnalyzer::UpdateFrame()
{
    int got, tocp,i,j,k,ii;
    int N = fsig_.N;
    int N2 = N/2;
    // int buflen = buflen;
    int analWinLen = fsig_.winsize/2;
    int synWinLen = analWinLen;
    float *ofp;                 /* RWD MUST be 32bit */
    float *fp;
    float *anal = analbuf;
    float *analOut = analbufOut;
    float *tempInput = input;
    float *analWindow = analwinbuf + analWinLen;
    float *tempOldInPhase = oldInPhase;
    float angleDif,real,imag,phase;
    float rratio;

    got = fsig_.overlap;      /*always assume */
    fp = overlapbuf;
    tocp = (got<= tempInput + buflen - nextIn ? got : tempInput + buflen - nextIn);
    got -= tocp;
    while (tocp-- > 0)
      *(nextIn++) = *fp++;

    if (got > 0) {
      nextIn -= buflen;
      while (got-- > 0)
        *nextIn++ = *fp++;
    }
    if (nextIn >= (tempInput + buflen))
      nextIn -= buflen;

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
    memset(anal, 0, sizeof(float)*(N+2));

    j = (nI - analWinLen - 1 + buflen) % buflen;     /*input pntr*/

    k = nI - analWinLen - 1;                 /*time shift*/
    while (k < 0)
      k += N;
    k = k % N;
    for (i = -analWinLen; i <= analWinLen; i++) {
      if (++j >= buflen)
        j -= buflen;
      if (++k >= N)
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
    
    fft_.Direct(anal, analOut, N);
    Interlace(analOut, N);

    /* conversion: The real and imaginary values in anal are converted to
       magnitude and angle-difference-per-second (assuming an
       intermediate sampling rate of rIn) and are returned in
       anal. */

    /*if (format==PVS_AMP_FREQ) {*/
    for (i=ii=0    /*,i0=anal,i1=anal+1,oi=tempOldInPhase*/;
         i <= N2;
         i++,ii+=2 /*i0+=2,i1+=2, oi++*/) {
      real = analOut[ii] /* *i0 */;
      imag = analOut[ii+1] /* *i1 */;
      /**i0*/ analOut[ii] = hypotf(real, imag);
      /* phase unwrapping */
      /*if (*i0 == 0.)*/
      if (/* *i0 */ analOut[ii] < 1.0E-10f)
        angleDif = 0.0f;
      else {
        rratio =  atan2((float)imag,(float)real);
        angleDif  = (phase = (float)rratio) - /**oi*/ tempOldInPhase[i];
        /* *oi */ tempOldInPhase[i] = phase;
      }

      if (angleDif > PI_F)
        angleDif = angleDif - TWOPI_F;
      if (angleDif < -PI_F)
        angleDif = angleDif + TWOPI_F;

      /* add in filter center freq.*/
      /* *i1 */ analOut[ii+1]  = angleDif * RoverTwoPi + ((float) i * Fexact);

    }
    /* } */
    /* else must be PVOC_COMPLEX */
    fp = analOut;
    ofp = fsig_.frame;      /* RWD MUST be 32bit */
    for (i=0;i < N+2;i++)
      /* *ofp++ = (float)(*fp++); */
      ofp[i] = (float) fp[i];

    nI += fsig_.overlap;                          /* increment time */
    if (nI > (synWinLen + fsig_.overlap))
      Ii = /*I*/fsig_.overlap;
    else
      if (nI > synWinLen)
        Ii = nI - synWinLen;
      else {
        Ii = 0;

        /*  for (i=nO+synWinLen; i<buflen; i++)
            if (i > 0)
            *(output+i) = 0.0f;
            */
      }

    IOi = Ii;
}

void SpectralAnalyzer::Interlace(float* fftSeparated, const int length) {
  // unfortunately, interleaving in place is not trivial, so another buffer will have to do
  int halflen = length / 2;
  for (int i = 0; i < halflen; i++) {
    swapBuffer_[i*2] = fftSeparated[i];
    swapBuffer_[i*2 + 1] = fftSeparated[i + halflen];
  }

  for (int i = 0; i < length; i++) {
    fftSeparated[i] = swapBuffer_[i];
  }

}