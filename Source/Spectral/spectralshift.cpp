#include <cmath>
#include <cstring>
#include "spectralshift.h"

using namespace daicsp;

void SpectralShift::Init(SpectralBuffer& fsig_in,
              float           shift,
              float           lowest,
              int             sample_rate,
              FORMANT         keepform,
              float           gain,
              int             coefficients)
{
  status_ = STATUS::OK;

  shift_ = shift;
  lowest_ = lowest;
  sample_rate_ = sample_rate;
  keepform_ = keepform;
  gain_ = gain;
  coefficients_ = coefficients;

  int32_t    N = fsig_in.N;

  if(fsig_in == fsig_out_)
  {
      status_ = STATUS::E_FSIG_EQUAL;
      return;
  }

  // if (UNLIKELY(p->fin == p->fout))
  //   csound->Warning(csound, Str("Unsafe to have same fsig as in and out"));
  if (fsig_in.sliding) {
    // if (fsig_out_.frame.auxp==NULL ||
    //     CS_KSMPS*(N+2)*sizeof(float) > (uint32_t)fsig_out_.frame.size)
    //   csound->AuxAlloc(csound, CS_KSMPS*(N+2)*sizeof(float),&fsig_out_.frame);
    // else memset(fsig_out_.frame.auxp, 0, CS_KSMPS*(N+2)*sizeof(float));

    status_ = STATUS::E_SLIDING_NOT_IMPLEMENTED;
    return;
  }
  // else
  //   {
  //     if (fsig_out_.frame.auxp == NULL ||
  //         fsig_out_.frame.size < sizeof(float) * (N + 2))  /* RWD MUST be 32bit */
  //       csound->AuxAlloc(csound, (N + 2) * sizeof(float), &fsig_out_.frame);
  //     else memset(fsig_out_.frame.auxp, 0, (N+2)*sizeof(float));
  //   }

  fsig_out_.N = N;
  fsig_out_.overlap = fsig_in.overlap;
  fsig_out_.winsize = fsig_in.winsize;
  fsig_out_.wintype = fsig_in.wintype;
  fsig_out_.format = fsig_in.format;
  fsig_out_.framecount = 1;
  lastframe_ = 0;
  fsig_out_.sliding = fsig_in.sliding;
  fsig_out_.NB = fsig_in.NB;
  fsig_out_.ready = false;

  // All equal to kFFTMaxFloats
  // if (p->ceps.auxp == NULL ||
  //     p->ceps.size < sizeof(float) * (N+2))
  //   csound->AuxAlloc(csound, sizeof(float) * (N + 2), &p->ceps);
  // else
  //   memset(p->ceps.auxp, 0, sizeof(float)*(N+2));
  // if (p->fenv.auxp == NULL ||
  //     p->fenv.size < sizeof(float) * (N+2))
  //   csound->AuxAlloc(csound, sizeof(float) * (N + 2), &p->fenv);
  // else
  //   memset(p->fenv.auxp, 0, sizeof(float)*(N+2));
  // if (p->ftmp.auxp == NULL ||
  //     p->ftmp.size < sizeof(float) * (N+4))
  //   csound->AuxAlloc(csound, sizeof(float) * (N + 2), &p->ftmp);

  // return OK;
}

SpectralBuffer& SpectralShift::ParallelProcess(SpectralBuffer& fsig_in)
{
  if (fsig_in.ready)
  {
    int32_t     i, chan, newchan, N = fsig_out_.N;
    float       pshift = shift_;
    int32_t     lowest = abs((int32_t) (lowest_ * N * (1.0f / sample_rate_)));
    int32_t     cshift = (int32_t) (pshift * N * (1.0f / sample_rate_));
    int32_t     keepform = (int32_t) keepform_;
    float   g = gain_;
    float   *fin = (float *) fsig_in.frame;
    float   *fout = (float *) fsig_out_.frame;
    float   *ftmp = (float *) ftmp_;
    float   *fenv = (float *) fenv_;
    // NOTE -- unused until keepform is implemented
    // float       max = 0.0f;
    // float   *ceps = (float *) ceps_;
    // float   binf;
    // int32_t coefs = (int32_t) coefficients_;

    // if (UNLIKELY(fout == NULL)) goto err1;
    // NOTE -- sliding not yet supported
    // if (fsig_in.sliding) {
    //   uint32_t offset = p->h.insdshead->ksmps_offset;
    //   uint32_t n, nsmps = CS_KSMPS;
    //   int32_t NB  = fsig_out_.NB;
    //   float g = *p->gain;
    //   lowest = lowest ? (lowest > NB ? NB : lowest) : 1;

    //   for (n=0; n<offset; n++) {
    //     CMPLX *fout = (CMPLX *) fsig_out_.frame.auxp + n*NB;
    //     for (i = 0; i < NB; i++) fout[i].re = fout[i].im = FL(0.0);
    //   }
    //   for (n=offset; n<nsmps; n++) {  float max = FL(0.0);
    //     CMPLX *fin = (CMPLX *) fsig_in.frame.auxp + n*NB;
    //     CMPLX *fout = (CMPLX *) fsig_out_.frame.auxp + n*NB;
    //     fout[0] = fin[0];
    //     fout[NB-1] = fin[NB-1];
    //     if (IS_ASIG_ARG(p->kshift)) {
    //       pshift = (p->kshift)[n];
    //     }
    //     for (i = 1; i < NB-1; i++) {
    //       if (keepform && (max < fin[i].re)) max = fin[i].re;
    //       if (i < lowest) {
    //         fout[i] = fin[i];
    //       }
    //     }
    //     for (i = lowest; i < NB; i++) {
    //       if (keepform == 0 || keepform == 1 || !max)
    //         fout[i].re = fin[i].re;
    //       else
    //         fout[i].re = fin[i].re * (fin[i].re / max);
    //       fout[i].im = (fin[i].im + pshift);
    //       /* Remove aliases */
    //       if (fout[i].im>=CS_ESR*0.5 ||
    //           fout[i].im<= -CS_ESR*0.5)
    //         fout[i].re = 0.0;
    //     }
    //     if (g!=1.0f)
    //       for (i = lowest; i < NB; i++) {
    //         fout[i].re *= g;
    //       }
    //   }
    //   return OK;
    // }
    if (lastframe_ < fsig_in.framecount) {
      int32_t j;
      lowest = lowest ? (lowest > N / 2 ? N / 2 : lowest << 1) : 2;

      fout[0] = fin[0];
      fout[N] = fin[N];
      memcpy(ftmp, fin, sizeof(float)*(N+2));

      for (j = i = 2; i < N; i += 2, j++) {
        fenv[j] = 0.0;
        if (i < lowest) {
          fout[i] = fin[i];
          fout[i + 1] = fin[i + 1];
        }
        else {
          fout[i] = 0.0f;
          fout[i + 1] = -1.0f;
        }
      }
      // NOTE -- keepform not supported
      if (keepform) {
        status_ = STATUS::W_KEEPFORM_NOT_IMPLEMENTED;
      }
      // if (keepform) { /* new modes 1 & 2 */
      //   int32_t cond = 1;
      //   int32_t tmp = N/2;
      //   tmp = tmp + tmp%2;
      //   for (i=j=0; i < N; i+=2, j++)
      //     fenv[j] = LOG(fin[i] > FL(0.0) ? fin[i] : FL(1e-20));
      //   if (coefs < 1) coefs = 80;
      //   while(cond) {
      //     cond = 0;
      //     for (j=i=0; i < N; i+=2, j++) {
      //       ceps[i] = fenv[j];
      //       ceps[i+1] = FL(0.0);
      //     }
      //     if (!(N & (N - 1)))
      //       csound->InverseComplexFFT(csound, ceps, N/2);
      //     else
      //       csoundInverseComplexFFTnp2(csound, ceps, tmp);
      //     for (i=coefs; i < N-coefs; i++) ceps[i] = 0.0;
      //     if (!(N & (N - 1)))
      //       csound->ComplexFFT(csound, ceps, N/2);
      //     else
      //       csoundComplexFFTnp2(csound, ceps, tmp);
      //     for (i=j=0; i < N; i+=2, j++) {
      //       if (keepform > 1) {
      //         if (fenv[j] < ceps[i])
      //           fenv[j] = ceps[i];
      //         if ((LOG(fin[i]) - ceps[i]) > 0.23) cond = 1;
      //       }
      //       else
      //         {
      //           fenv[j] = EXP(ceps[i]);
      //           max = max < fenv[j] ? fenv[j] : max;
      //         }
      //     }
      //   }
      //   if (keepform > 1)
      //     for (i=0; i<N/2; i++) {
      //       fenv[i] = EXP(fenv[i]);
      //       max = max < fenv[i] ? fenv[i] : max;
      //     }
      //   if (max)
      //     for (j=i=lowest; i<N; i+=2, j++) {
      //       fenv[j]/=max;
      //       binf = (j)*sample_rate_/N;
      //       if (fenv[j] && binf < sample_rate_/2+pshift )
      //         ftmp[i] /= fenv[j];
      //     }
      // }

      // if(keepform) 
      // {
      //   for (i = lowest, chan = lowest >> 1; i < N; chan++, i += 2) {
      //     newchan = (chan + cshift) << 1;
      //     if (newchan < N && newchan > lowest) {
      //       fout[newchan] = ftmp[i] * fenv[newchan>>1];
      //       fout[newchan + 1] = (float) (ftmp[i + 1] + pshift);
      //     }
      //   }
      // } else 
      // {
        for (i = lowest, chan = lowest >> 1; i < N; chan++, i += 2) {
          newchan = (chan + cshift) << 1;
          if (newchan < N && newchan > lowest) {
            fout[newchan] = ftmp[i];
            fout[newchan + 1] = (float) (ftmp[i + 1] + pshift);
          }
        }

      // }

      for (i = lowest; i < N; i += 2) {
        if (fout[i + 1] == -1.0f)
          fout[i] = 0.0f;
        else
          fout[i] *= g;
      }

      fsig_out_.framecount = lastframe_ = fsig_in.framecount;
    }
  }

  fsig_out_.ready = fsig_in.ready;
  return fsig_out_;
//   return OK;
//  err1:
//   return csound->PerfError(csound, &(p->h),
//                            Str("pvshift: not initialised"));
}