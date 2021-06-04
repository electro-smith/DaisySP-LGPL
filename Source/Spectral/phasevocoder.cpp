#include <math.h>
#include <string.h>
#include "dsp.h"
#include "phasevocoder.h"

using namespace daicsp;

void PhaseVocoder::Init(float sampleRate)
{
    // Set private members to defaults
    sr_ = sampleRate;
    // Do stuff
    
}

void PhaseVocoder::Process(const float &in)
{
    int offset = h.insdshead->ksmps_offset;
    int early  = h.insdshead->ksmps_no_end;
    int i, nsmps = CS_KSMPS;
    float *aout = aout;

    // if (output.auxp==NULL) {
    //   return csound->PerfError(csound,&(h),
    //                            Str("pvsynth: Not Initialised.\n"));
    // }
    if (fsig->sliding) return Analyze();
    if (offset) memset(aout, '\0', offset*sizeof(float));
    if (early) {
      nsmps -= early;
      memset(&aout[nsmps], '\0', early*sizeof(float));
    }
    for (i=offset; i<nsmps; i++)
      aout[i] = Tick();
    // return OK;
}

void PhaseVocoder::Analyze(float sample) 
{
    int i, k;
    int ksmps = CS_KSMPS;
    int N = fsig->N;
    int NB = fsig->NB;
    float *aout = aout;
    Complex *ff;
    float *h = (float*)oldOutPhase.auxp;
    float *output = (float*)output.auxp;

    /* Get real part from AMP/FREQ */
    for (i=0; i<ksmps; i++) {
      float a;
      ff = (Complex*)(fsig->frame.auxp) + i*NB;
      for (k=0; k<NB; k++) {
        float tmp, phase;

        tmp = ff[k].b; /* Actually frequency */
        /* subtract bin mid frequency */
        tmp -= (float)k * csound->esr/N;
        /* get bin deviation from freq deviation */
        tmp *= TWOPI_F /csound->esr;
        /* add the overlap phase advance back in */
        tmp += (float)k*TWOPI_F/N;
        h[k] = phase = mod2Pi(h[k] + tmp);
        output[k] = ff[k].a*cos(phase);
      }
      a = 0.0f;
      for (k=1; k<NB-1; k++) {
        a -= output[k];
        if (k+1<NB-1) a+=output[++k];
      }
      aout[i] = (a+a+output[0]-output[NB-1])/N;
    }
    // return OK;
}

void PhaseVocoder::Tick(float sample) 
{
    float *outbuf = (float *) (overlapbuf.auxp);

    if (outptr== fsig->overlap) {
      UpdateFrame();
      outptr = 0;
    }
    return outbuf[outptr++];
}

void PhaseVocoder::UpdateFrame()
{
    int i,j,k,ii,NO,NO2;
    float *anal;                                        /* RWD MUST be 32bit */
    float *syn, *output;
    float *oldOutPhase = (float *) (oldOutPhase.auxp);
    int N = fsig->N;
    float *obufptr,*outbuf,*synWindow;
    float mag,phase,angledif, the_phase;
    int synWinLen = fsig->winsize / 2;
    int overlap = fsig->overlap;
    /*int32 format = fsig->format; */

    /* fsigs MUST be corect format, as we offer no mechanism for
       assignment to a different one*/

    NO = N;        /* always the same */
    NO2 = NO/2;
    syn = (float *) (synbuf.auxp);
    anal = (float *) (fsig->frame.auxp);             /* RWD MUST be 32bit */
    output = (float *) (output.auxp);
    outbuf = (float *) (overlapbuf.auxp);
    synWindow = (float *) (synwinbuf.auxp) + synWinLen;

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
    the_phase = /* *(oldOutPhase + i) */ oldOutPhase[i] + angledif;
    if (i== bin_index)
        the_phase = (float) fmod(the_phase,TWOPI_F);
    /* *(oldOutPhase + i) = the_phase; */
    oldOutPhase[i] = the_phase;
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
    if (!(NO & (NO - 1))) {
      /*printf("N %d %d \n", NO, NO & (NO-1));*/
      syn[1] = syn[NO];
      /* csound->InverseRealFFT(csound, syn, NO);*/
      csound->RealFFT2(csound,setup,syn);
      syn[NO] = syn[NO + 1] = 0.0f;
    }
    else
      csound->InverseRealFFTnp2(csound, syn, NO);
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
