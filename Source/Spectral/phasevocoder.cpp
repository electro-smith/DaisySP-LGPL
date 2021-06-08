#include <math.h>
#include <string.h>
#include "dsp.h"
#include "phasevocoder.h"
#include "spectral.h"

using namespace daicsp;

void PhaseVocoder::Init(SpectralBuffer fsig, size_t sampleRate, size_t block)
{
    sr_ = sampleRate;

    float *analwinhalf;
    float *synwinhalf;
    float sum;
    int halfwinsize;
    // int i,nBins,Mf,Lf;
    int i, Mf, Lf;
    double IO;

    /* get params from input fsig */
    /* we TRUST they are legal */
    int N = fsig.N;
    int overlap = fsig.overlap;
    int M = fsig.winsize;
    int wintype = fsig.wintype;

    // TODO -- ensure these are not necessary
    // p->fftsize = N;
    // p->winsize = M;
    // p->overlap = overlap;
    // p->wintype = wintype;
    // p->format = fsig.format;

    // NOTE -- Sliding is always set to false, so this will never happen
    // if (fsig.sliding) {
    //   /* get params from input fsig */
    //   /* we TRUST they are legal */
    //   int wintype = fsig.wintype;
    //   /* and put into locals */
    //   p->wintype = wintype;
    //   p->format = fsig.format;
    //   csound->AuxAlloc(csound, fsig.NB * sizeof(double), &p->oldOutPhase);
    //   csound->AuxAlloc(csound, fsig.NB * sizeof(double), &p->output);
    //   return OK;
    // }

    /* and put into locals */
    halfwinsize = M/2;
    buflen = M*4;
    IO = (double)overlap;         /* always, no time-scaling possible */

    // TODO -- check what csound->esr is
    // float arate = csound->esr / (float) overlap;
    // float fund = csound->esr / (float) N;

    float arate = sr_ / (float) overlap;
    // float fund = sr_ / (float) N;

    // nBins = N/2 + 1;

    Lf = Mf = 1 - M%2;

    // TODO -- take note of these and how they may affect static allocation
    // /* deal with iinit later on! */
    // csound->AuxAlloc(csound, overlap * sizeof(float), &p->overlapbuf);
    // csound->AuxAlloc(csound, (N+2) * sizeof(float), &p->synbuf);
    // csound->AuxAlloc(csound, (M+Mf) * sizeof(float), &p->analwinbuf);
    // csound->AuxAlloc(csound, (M+Mf) * sizeof(float), &p->synwinbuf);
    // csound->AuxAlloc(csound, nBins * sizeof(float), &p->oldOutPhase);
    // csound->AuxAlloc(csound, buflen * sizeof(float), &p->output);

    synwinhalf = synwinbuf + halfwinsize;

    /* synthesis windows */
    if (M <= N) {
      SpectralWindow(synwinhalf, wintype, M);

      for (i = 1; i <= halfwinsize; i++)
        *(synwinhalf - i) = *(synwinhalf + i - Lf);


       sum = 0.0f;
        for (i = -halfwinsize; i <= halfwinsize; i++)
           sum += *(synwinhalf + i);
         sum = 2.0f / sum;

       for (i = -halfwinsize; i <= halfwinsize; i++)
        *(synwinhalf + i) *= sum;

         sum = 0.0f;
        /* no timescaling, so I(nterpolation) will always = D(ecimation) = overlap */
        for (i = -halfwinsize; i <= halfwinsize; i+=overlap)
          sum += *(synwinhalf + i) * *(synwinhalf + i);
    }
    else {
     /* have to make analysis window to get amp scaling */
    /* so this ~could~ be a local alloc and free...*/
      double dN = (double)N;
      analwinhalf = analwinbuf + halfwinsize;
    SpectralWindow(analwinhalf, wintype, M);

    for (i = 1; i <= halfwinsize; i++){
      analwinhalf[-i] = analwinhalf[i - Mf];
    }

      // sinc function
    if (Mf) {
        *analwinhalf *= (float)(dN * sin(HALFPI_F/dN) / ( HALFPI_F));
    }
      for (i = 1; i <= halfwinsize; i++)
        *(analwinhalf + i) *= (float)
          (dN * sin((double)(PI_F*(i+0.5*Mf)/dN)) / (PI_F*(i+0.5*Mf)));
      for (i = 1; i <= halfwinsize; i++)
        *(analwinhalf - i) = *(analwinhalf + i - Mf);

     /* get net amp */
    sum = 0.0f;
    for (i = -halfwinsize; i <= halfwinsize; i++)
    sum += *(analwinhalf + i);
    sum = 2.0f / sum;  /* factor of 2 comes in later in trig identity */

      SpectralWindow(synwinhalf, wintype, M);

      for (i = 1; i <= halfwinsize; i++)
        *(synwinhalf - i) = *(synwinhalf + i - Lf);

      if (Lf)
        *synwinhalf *= (float)(IO * sin((double)(HALFPI_F/IO)) / (HALFPI_F));
      for (i = 1; i <= halfwinsize; i++)
        *(synwinhalf + i) *= (float)
          ((double)IO * sin((double)(PI_F*(i+0.5*Lf)/IO)) /
           (PI_F*(i+0.5*(double)Lf)));
      for (i = 1; i <= halfwinsize; i++)
        *(synwinhalf - i) = *(synwinhalf + i - Lf);
    }

    // unless casting structs to void does something magical,
    // csound->GetInverseRealFFTScale just returns 1
//    if (!(N & (N - 1L)))
//      sum = csound->GetInverseRealFFTScale(csound, (int) N)/ sum;
//     else
//       sum = 1.0f / sum;

    sum = 1.0f / sum;

  for (i = -halfwinsize; i <= halfwinsize; i++)
      *(synwinhalf + i) *= sum;

/*  p->invR = FL(1.0) / csound->esr; */
    RoverTwoPi = arate / TWOPI_F;
    TwoPioverR = TWOPI_F / arate;
    // Fexact =  csound->esr / (float)N;
    Fexact = sr_ / (float) N;
    nO = -(halfwinsize / overlap) * overlap; /* input time (in samples) */
    Ii = 0;                          /* number of new outputs to write */
    IOi = 0;
    outptr = 0;
    nextOut = output;

    // if (!(N & (N - 1))) /* if pow of two use this */
    //   p->setup = csound->RealFFT2Setup(csound,N,FFT_INV);
    // return OK;

    fft_.Init();
}

float* PhaseVocoder::Process(SpectralBuffer &fsig, size_t size)
{
    // int offset = h.insdshead->ksmps_offset;
    // int early  = h.insdshead->ksmps_no_end;
    int offset = 0;
    // int early = 0;

    int i, nsmps = size;

    // if (output.auxp==NULL) {
    //   return csound->PerfError(csound,&(h),
    //                            Str("pvsynth: Not Initialised.\n"));
    // }
    
    // // This will always be false
    // if (fsig.sliding) return Analyze();

    // TODO -- could these be useful for handling 
    // if (offset) memset(outputBuffer, '\0', offset*sizeof(float));
    // if (early) {
    //   nsmps -= early;
    //   memset(&outputBuffer[nsmps], '\0', early*sizeof(float));
    // }
    for (i=offset; i<nsmps; i++)
      outputBuffer[i] = Tick(fsig);
    // return OK;
    return outputBuffer;
}

void PhaseVocoder::Analyze(float sample) 
{
    // int i, k;
    // int ksmps = CS_KSMPS;
    // int N = fsig.N;
    // int NB = fsig.NB;
    // float *outputBuffer = outputBuffer;
    // Complex *ff;
    // float *h = (float*)oldOutPhase.auxp;
    // float *output = (float*)output.auxp;

    // /* Get real part from AMP/FREQ */
    // for (i=0; i<ksmps; i++) {
    //   float a;
    //   ff = (Complex*)(fsig.frame.auxp) + i*NB;
    //   for (k=0; k<NB; k++) {
    //     float tmp, phase;

    //     tmp = ff[k].b; /* Actually frequency */
    //     /* subtract bin mid frequency */
    //     tmp -= (float)k * csound->esr/N;
    //     /* get bin deviation from freq deviation */
    //     tmp *= TWOPI_F /csound->esr;
    //     /* add the overlap phase advance back in */
    //     tmp += (float)k*TWOPI_F/N;
    //     h[k] = phase = mod2Pi(h[k] + tmp);
    //     output[k] = ff[k].a*cos(phase);
    //   }
    //   a = 0.0f;
    //   for (k=1; k<NB-1; k++) {
    //     a -= output[k];
    //     if (k+1<NB-1) a+=output[++k];
    //   }
    //   outputBuffer[i] = (a+a+output[0]-output[NB-1])/N;
    // }
    // // return OK;
}

float PhaseVocoder::Tick(SpectralBuffer& fsig) 
{
    if (outptr == fsig.overlap) {
      UpdateFrame(fsig);
      outptr = 0;
    }
    return overlapbuf[outptr++];
}

void PhaseVocoder::UpdateFrame(SpectralBuffer& fsig)
{
    int i,j,k,ii,NO,NO2;
    float *anal;                                        /* RWD MUST be 32bit */
    float *syn;
    float *tempOldOutPhase = oldOutPhase;
    int N = fsig.N;
    float *obufptr,*outbuf,*synWindow;
    float mag,phase,angledif, the_phase;
    int synWinLen = fsig.winsize / 2;
    int overlap = fsig.overlap;
    /*int32 format = fsig.format; */

    /* fsigs MUST be corect format, as we offer no mechanism for
       assignment to a different one*/

    NO = N;        /* always the same */
    NO2 = NO/2;
    syn = synbuf;
    anal = fsig.frame;             /* RWD MUST be 32bit */
    // output = (float *) (output.auxp);
    outbuf = overlapbuf;
    synWindow = synwinbuf + synWinLen;

    /* reconversion: The magnitude and angle-difference-per-second in syn
       (assuming an intermediate sampling rate of rOut) are
       converted to real and imaginary values and are returned in syn.
       This automatically incorporates the proper phase scaling for
       time modifications. */

    if (NO <= N) {
      for (i = 0; i < NO+2; i++)
        syn[i] = (float) anal[i];
    }
    else {
      for (i = 0; i <= N+1; i++)
        syn[i] = (float) anal[i];
      for (i = N+2; i < NO+2; i++)
        syn[i] = 0.0f;
    }

    for (i=ii=0 /*, i0=syn, i1=syn+1*/; i<= NO2; i++, ii+=2 /*i0+=2,  i1+=2*/) {
    mag = syn[ii]; /* *i0; */
    /* RWD variation to keep phase wrapped within +- TWOPI_F */
    /* this is spread across several frame cycles, as the problem does not
        develop for a while */

    angledif = TwoPioverR * ( /* *i1 */ syn[ii+1] - ((float)i * Fexact));
    the_phase = /* *(oldOutPhase + i) */ tempOldOutPhase[i] + angledif;
    if (i== bin_index)
        the_phase = (float) fmod(the_phase,TWOPI_F);
    /* *(oldOutPhase + i) = the_phase; */
    tempOldOutPhase[i] = the_phase;
    phase = the_phase;
    /* *i0 */ syn[ii]  = (float)((float)mag * cos((float)phase));
    /* *i1 */ syn[ii+1] = (float)((float)mag * sin((float)phase));
    }

    /* for phase normalization */
    if (++(bin_index) == NO2+1)
      bin_index = 0;

    /* else it must be PVOC_COMPLEX */

    /* synthesis: The synthesis subroutine uses the Weighted Overlap-Add
       technique to reconstruct the time-domain signal.  The (N/2 + 1)
       phase vocoder channel outputs at time n are inverse Fourier
       transformed, windowed, and added into the output array.  The
       subroutine thinks of output as a shift register in which
       locations are referenced modulo obuflen.  Therefore, the main
       program must take care to zero each location which it "shifts"
       out (to standard output). The subroutines reals and fft
       together perform an efficient inverse FFT.  */
    // if (!(NO & (NO - 1))) {
    //   /*printf("N %d %d \n", NO, NO & (NO-1));*/
    //   syn[1] = syn[NO];
    //   /* csound->InverseRealFFT(csound, syn, NO);*/
    //   csound->RealFFT2(csound,setup,syn);
    //   syn[NO] = syn[NO + 1] = 0.0f;
    // }
    // else
    //   csound->InverseRealFFTnp2(csound, syn, NO);

    Deinterlace(syn, synbufOut, NO);
    if (NO != WINDOW_SIZE::MAX)
    {
      int num_passes = GetPasses(NO);
      fft_.Inverse(syn, synbufOut, num_passes);
    }
    else
    {
      fft_.Inverse(syn, synbufOut);
    }
    // NOTE -- shy_fft outputs raw inverse values, ranging from +- the 
    // size of the fft. Also, it outputs to a second buffer, but we
    // want to use the first one (syn)
    for (int n = 0; n < NO; n++) {
        syn[n] = synbufOut[n] / NO;
    }

    j = nO - synWinLen - 1;
    while (j < 0)
      j += buflen;
    j = j % buflen;

    k = nO - synWinLen - 1;
    while (k < 0)
      k += NO;
    k = k % NO;

    for (i = -synWinLen; i <= synWinLen; i++) { /*overlap-add*/
      if (++j >= buflen)
        j -= buflen;
      if (++k >= NO)
        k -= NO;
      /* *(output + j) += *(syn + k) * *(synWindow + i); */
      output[j] += syn[k] * synWindow[i];
    }

    obufptr = outbuf;

    for (i = 0; i < IOi;) {  /* shift out next IOi values */
      int todo = (IOi-i <= output+buflen - nextOut ?
                  IOi-i : output+buflen - nextOut);
      /*outfloats(nextOut, todo, ofd);*/
      /*copy data to external buffer */
      /*for (n=0;n < todo;n++)
       *obufptr++ = nextOut[n]; */
      memcpy(obufptr, nextOut, sizeof(float)*todo);
      obufptr += todo;

      i += todo;

      /* for (j = 0; j < todo; j++)
       *nextOut++ = 0.0f; */
      memset(nextOut, 0, sizeof(float)*todo);
      nextOut += todo;

      if (nextOut >= (output + buflen))
        nextOut -= buflen;
    }

    /* increment time */
    nO += overlap;

    if (nO > (synWinLen + /*I*/overlap))
      Ii = overlap;
    else
      if (nO > synWinLen)
        Ii = nO - synWinLen;
      else {
        Ii = 0;
        for (i=nO+synWinLen; i<buflen; i++)
          if (i > 0)
            output[i] = 0.0f;
      }
    IOi =  Ii;
}

void PhaseVocoder::Deinterlace(float* interlaced, float* workingBuffer, const int length) 
{
    int halflen = length / 2;

    for (int i = 0; i < halflen; i++) {
        workingBuffer[i] = interlaced[i * 2];
        workingBuffer[i + halflen] = interlaced[i * 2 + 1];
    }

    for (int i = 0; i < length; i++) {
        interlaced[i] = workingBuffer[i];
    }
}
