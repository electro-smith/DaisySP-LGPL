#include <cstddef>
#include <cstring>
#include <cmath>
#include "spectralscale.h"

using namespace daicsp;
using namespace std;

void SpectralScale::Init(SpectralBuffer& fsig_in,
                         float           scale,
                         int             sample_rate,
                         FORMANT         formants,
                         float           gain,
                         int             coefficients)
{
    status_      = STATUS::OK;
    sample_rate_ = sample_rate;

    kscal_    = scale;
    keepform_ = formants;
    gain_     = gain;
    coefs_    = coefficients;

    int N = fsig_in.N, tmp;

    //   if (UNLIKELY(fin == fsig_out_))
    //     csound->Warning(csound, Str("Unsafe to have same fsig as in and out"));
    if(fsig_in == fsig_out_)
    {
        status_ = STATUS::E_FSIG_EQUAL;
        return;
    }
    fsig_out_.NB      = fsig_in.NB;
    fsig_out_.sliding = fsig_in.sliding;
    if(fsig_in.sliding)
    {
        status_ = STATUS::E_SLIDING_NOT_IMPLEMENTED;
        return;

        // if (fsig_out_.frame.auxp == NULL ||
        //     fsig_out_.frame.size < CS_KSMPS * sizeof(float) * (N + 2))
        //   csound->AuxAlloc(csound, CS_KSMPS * sizeof(float) * (N + 2),
        //                    &fsig_out_.frame);
    }
    else
    {
        //   if (fsig_out_.frame.auxp == NULL ||
        //       fsig_out_.frame.size < sizeof(float) * (N + 2))  /* RWD MUST be 32bit */
        //     csound->AuxAlloc(csound, sizeof(float) * (N + 2), &fsig_out_.frame);
    }

    //   if (ftmp.auxp == NULL ||
    //       ftmp.size < sizeof(float) * (N+4))
    //     csound->AuxAlloc(csound, sizeof(float) * (N + 2), &ftmp);

    fsig_out_.N          = N;
    fsig_out_.overlap    = fsig_in.overlap;
    fsig_out_.winsize    = fsig_in.winsize;
    fsig_out_.wintype    = fsig_in.wintype;
    fsig_out_.format     = fsig_in.format;
    fsig_out_.framecount = 1;
    fsig_out_.ready      = false;
    lastframe            = 0;
    tmp                  = N + N % 2;
    //   if (ceps.auxp == NULL ||
    //       ceps.size < sizeof(float) * (tmp+2))
    //     csound->AuxAlloc(csound, sizeof(float) * (tmp + 2), &ceps);
    memset(ceps_, 0, sizeof(float) * (tmp + 2));
    //   if (fenv.auxp == NULL ||
    //       fenv.size < sizeof(float) * (N+2))
    //     csound->AuxAlloc(csound, sizeof(float) * (N + 2), &fenv);
    memset(fenv_, 0, sizeof(float) * (N + 2));
    //   fwdsetup = csound->RealFFT2Setup(csound, N/2, FFT_FWD);
    //   invsetup = csound->RealFFT2Setup(csound, N/2, FFT_INV);
    //   return OK;
    fft_.Init();
}

SpectralBuffer& SpectralScale::ParallelProcess(SpectralBuffer& fsig_in)
{
    if(fsig_in.ready)
    {
        int32_t i, chan, N = fsig_out_.N;
        float   max = 0.0f;
        //   float   *fin = fsig_in.frame;
        //   float   *fout = fsig_out_.frame;

        // ensuring audio callback doesn't fire during assignment
        float   pscal;
        int32_t keepform;
        float   g;
        int32_t coefs;
        {
            // TODO -- add a blocker to prevent parameter mis-match
            // daisy::ScopedIrqBlocker block;

            pscal    = fabs(kscal_);
            keepform = keepform_;
            g        = gain_;
            coefs    = coefs_;
        }

        float* fenv = fenv_;
        float* ftmp = ftmp_;
        float* ceps = ceps_;

        float sr = sample_rate_, binf;

        // NOTE -- sliding not yet implemented
        //   if (fsig_out_.sliding) {
        //     uint32_t offset = h.insdshead->ksmps_offset;
        //     uint32_t n, nsmps = CS_KSMPS;
        //     int32_t NB    = fsig_out_.NB;
        //     float   g = *gain;
        //     for (n=0; n<offset; n++) {
        //       CMPLX   *fsig_out_ = (CMPLX *) fsig_out_.frame.auxp + n*NB;
        //       for (i = 0; i < NB; i++) fsig_out_[i].re = fsig_out_[i].im = FL(0.0);
        //     }
        //     for (n=offset; n<nsmps; n++) {
        //       float    max = FL(0.0);
        //       CMPLX   *fin = (CMPLX *) fin->frame.auxp + n*NB;
        //       CMPLX   *fsig_out_ = (CMPLX *) fsig_out_.frame.auxp + n*NB;

        //       fsig_out_[0] = fin[0];
        //       fsig_out_[NB-1] = fin[NB-1];
        //       if (IS_ASIG_ARG(kscal)) {
        //         pscal = FABS(kscal[n]);
        //       }
        //       if (keepform)
        //         for (i = 1; i < NB-1; i++) {
        //           max = max < fin[i].re ? fin[i].re : max;
        //         }

        //       for (i = 1; i < NB-1; i++) {
        //         if (keepform == 0 || keepform == 1 || !max)
        //           fsig_out_[i].re = fin[i].re;
        //         else
        //           fsig_out_[i].re = fin[i].re * (fin[i].re / max);
        //         fsig_out_[i].im = fin[i].im * pscal;
        //         /* Remove aliases */
        //         if (fsig_out_[i].im>=CS_ESR*0.5 ||
        //             fsig_out_[i].im<= -CS_ESR*0.5)
        //           fsig_out_[i].re=0.0;
        //       }

        //       for (i = 1; i < NB; i++) {
        //         fsig_out_[i].re *= g;
        //       }
        //     }
        //     return OK;
        //   }

        if(lastframe < fsig_in.framecount)
        {
            int32_t n;
            fsig_out_.frame[0] = fsig_in.frame[0];
            fsig_out_.frame[N] = fsig_in.frame[N];
            memcpy(ftmp, fsig_in.frame, sizeof(float) * (N + 2));

            for(i = 2, n = 1; i < N; i += 2, n++)
            {
                fsig_out_.frame[i]     = 0.0f;
                fsig_out_.frame[i + 1] = -1.0f;
                fenv[n]                = 0.f;
            }

            if(keepform)
            {
                int32_t cond = 1;
                int32_t j;
                for(i = j = 0; i < N; i += 2, j++)
                    fenv[j] = log(ftmp[i] > 0.0f ? ftmp[i] : 1e-20f);


                if(keepform > 2)
                { /* experimental mode 3 */
                    int32_t w = 5, w2 = w * 2;
                    for(i = 0; i < w; i++)
                        ceps[i] = fenv[i];
                    for(i = w; i < N / 2 - w; i++)
                    {
                        ceps[i] = 0.0;
                        for(j = -w; j < w; j++)
                            ceps[i] += fenv[i + j];
                        ceps[i] /= w2;
                    }
                    for(i = 0; i < N / 2; i++)
                    {
                        fenv[i] = exp(ceps[i]);
                        max     = max < fenv[i] ? fenv[i] : max;
                    }
                    if(max)
                        for(j = i = 0; i < N; i += 2, j++)
                        {
                            fenv[j] /= max;
                            binf = (j)*sr / N;
                            if(fenv[j] && binf < pscal * sr / 2)
                                ftmp[i] /= fenv[j];
                        }
                }
                else
                { /* new modes 1 & 2 */
                    int32_t tmp = N / 2, j;
                    tmp         = tmp + tmp % 2;
                    if(coefs_ < 1)
                        coefs_ = 80;
                    while(cond)
                    {
                        cond = 0;
                        for(i = 0; i < N / 2; i++)
                        {
                            ceps[i] = fenv[i];
                        }

                        //   if (!(N & (N - 1))) // -> if N is a power of 2
                        //     csound->RealFFT2(csound, fwdsetup, ceps);
                        //   else
                        //     csound->RealFFTnp2(csound, ceps, tmp);

                        if(N != kFFTMaxSize)
                        {
                            int num_passes = GetPasses(tmp);
                            fft_.Direct(ceps, ceps_out_, num_passes);
                        }
                        else
                        {
                            fft_.Direct(ceps, ceps_out_);
                        }
                        Interleave(ceps_out_, ceps, N);

                        for(i = coefs; i < N / 2; i++)
                            ceps[i] = 0.0;

                        //   // NOTE -- is this a mistake?
                        //   if (!(N & (N - 1)))
                        //     csound->RealFFT2(csound, invsetup, ceps);
                        //   else
                        //     csound->InverseRealFFTnp2(csound, ceps, tmp);

                        Deinterleave(ceps, ceps_out_, N);
                        if(N != kFFTMaxSize)
                        {
                            int num_passes = GetPasses(N);
                            fft_.Inverse(ceps, ceps_out_, num_passes);
                        }
                        else
                        {
                            fft_.Inverse(ceps, ceps_out_);
                        }
                        for(int n = 0; n < N; n++)
                        {
                            ceps[n] = ceps_out_[n] / N;
                        }

                        for(i = j = 0; i < N / 2; i++, j += 2)
                        {
                            if(keepform > 1)
                            {
                                if(fenv[i] < ceps[i])
                                    fenv[i] = ceps[i];
                                if((log(ftmp[j]) - ceps[i]) > 0.23f)
                                    cond = 1;
                            }
                            else
                            {
                                fenv[i] = exp(ceps[i]);
                                max     = max < fenv[i] ? fenv[i] : max;
                            }
                        }
                    }
                    if(keepform > 1)
                        for(i = 0; i < N / 2; i++)
                        {
                            fenv[i] = exp(ceps[i]);
                            max     = max < fenv[i] ? fenv[i] : max;
                        }

                    if(max)
                        for(i = j = 2; i < N / 2; i++, j += 2)
                        {
                            fenv[i] /= max;
                            binf = (i)*sr / N;
                            if(fenv[i] && binf < pscal * sr / 2)
                                ftmp[j] /= fenv[i];
                        }
                }
            }
            if(keepform)
            {
                for(i = 2, chan = 1; i < N; chan++, i += 2)
                {
                    int32_t newchan;
                    newchan = (int32_t)((chan * pscal) + 0.5) << 1;
                    if(newchan < N && newchan > 0)
                    {
                        fsig_out_.frame[newchan] = ftmp[i] * fenv[newchan >> 1];
                        fsig_out_.frame[newchan + 1]
                            = (float)(ftmp[i + 1] * pscal);
                    }
                }
            }
            else
            {
                for(i = 2, chan = 1; i < N; chan++, i += 2)
                {
                    int32_t newchan;
                    newchan = (int32_t)((chan * pscal) + 0.5) << 1;
                    if(newchan < N && newchan > 0)
                    {
                        fsig_out_.frame[newchan] = ftmp[i];
                        fsig_out_.frame[newchan + 1]
                            = (float)(ftmp[i + 1] * pscal);
                    }
                }
            }

            for(i = 2; i < N; i += 2)
            {
                if(isnan(fsig_out_.frame[i]))
                    fsig_out_.frame[i] = 0.0f;
                if(fsig_out_.frame[i + 1] == -1.0f)
                {
                    fsig_out_.frame[i] = 0.f;
                }
                else
                    fsig_out_.frame[i] *= g;
            }
            fsig_out_.framecount = lastframe = fsig_in.framecount;
        }
    }

    fsig_out_.ready = fsig_in.ready;
    return fsig_out_;
    //   return OK;
}