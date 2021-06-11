#include "spectral.h"
#include "spectralblur.h"

using namespace daicsp;

void SpectralBlur::Init(SpectralBuffer<FFT::MAX_FLOATS>& fsigIn,
                        float                            delay,
                        float*                           delayBuffer,
                        size_t                           delaySize,
                        int                              sampleRate,
                        int                              block)
{
    status_   = STATUS::OK;
    sr_       = sampleRate;
    kdel_     = delay;
    delframes = delayBuffer;

    float* delayPtr;
    int    N    = fsigIn.N, i, j;
    int    olap = fsigIn.overlap;
    int    delayframes, framesize = N + 2;

    // Calculating the max delay time from the given buffer
    maxdel = (float)delaySize / ((N + 2) * ((float)sampleRate / olap));

    if(maxdel < delay)
    {
        status_ = STATUS::E_BUFFER_TOO_SMALL;
        return;
    }

    //   if (fsigIn == fsigOut_)
    //     csound->Warning(csound, Str("Unsafe to have same fsig as in and out"));
    if(&fsigIn == &fsigOut_)
    {
        status_ = STATUS::E_FSIG_EQUAL;
        return;
    }

    if(fsigIn.sliding)
    {
        status_ = STATUS::E_SLIDING_NOT_IMPLEMENTED;
        return;
        // csound->InitError(csound, Str("pvsblur does not work sliding yet"));
        // delayframes = (int) (FL(0.5) + *maxdel * CS_ESR);
        // if (fsigOut_.frame.auxp == NULL ||
        //     fsigOut_.frame.size < sizeof(MYFLT) * CS_KSMPS * (N + 2))
        //   csound->AuxAlloc(csound, (N + 2) * sizeof(MYFLT) * CS_KSMPS,
        //                    &fsigOut_.frame);

        // if (delframes.auxp == NULL ||
        //     delframes.size < (N + 2) * sizeof(MYFLT) * CS_KSMPS * delayframes)
        //   csound->AuxAlloc(csound,
        //                    (N + 2) * sizeof(MYFLT) * CS_KSMPS * delayframes,
        //                    &delframes);
    }
    else
    {
        frpsec      = (float)sr_ / olap;
        delayframes = (int)(maxdel * frpsec);

        //   NOTE -- No need for allocation with this scheme
        //   if (fsigOut_.frame.auxp == NULL ||
        //       fsigOut_.frame.size < sizeof(float) * (N + 2))
        //     csound->AuxAlloc(csound, (N + 2) * sizeof(float), &fsigOut_.frame);

        //   if (delframes.auxp == NULL ||
        //       delframes.size < (N + 2) * sizeof(float) * CS_KSMPS * delayframes)
        //     csound->AuxAlloc(csound, (N + 2) * sizeof(float) * delayframes,
        //                      &delframes);
    }
    delayPtr = delframes;

    for(j = 0; j < framesize * delayframes; j += framesize)
        for(i = 0; i < N + 2; i += 2)
        {
            delayPtr[i + j]     = 0.0f;
            delayPtr[i + j + 1] = i * (float)sr_ / N;
        }

    fsigOut_.N          = N;
    fsigOut_.overlap    = olap;
    fsigOut_.winsize    = fsigIn.winsize;
    fsigOut_.wintype    = fsigIn.wintype;
    fsigOut_.format     = fsigIn.format;
    fsigOut_.framecount = 1;
    lastframe           = 0;
    count               = 0;
    fsigOut_.sliding    = fsigIn.sliding;
    fsigOut_.NB         = fsigIn.NB;
}

SpectralBuffer<FFT::MAX_FLOATS>&
SpectralBlur::Process(SpectralBuffer<FFT::MAX_FLOATS>& fsigIn, int block)
{
    int    j, i, N = fsigOut_.N, first, framesize = N + 2;
    int    countr = count;
    float  amp = 0.0, freq = 0.0;
    int    delayframes = (int)(kdel_ * frpsec);
    int    kdel        = delayframes * framesize;
    int    mdel        = (int)(maxdel * frpsec) * framesize;
    float* fin         = fsigIn.frame;
    float* fout        = fsigOut_.frame;
    float* delay       = delframes;

    // if (UNLIKELY(fsigOut_ == NULL || delay == NULL)) goto err1;

    // NOTE -- sliding doesn't work yet
    // if (fsigIn.sliding) {
    //     unsigned int offset = 0;
    //     unsigned int n, nsmps = block;
    //     int NB = fsigIn.NB;
    //     kdel = kdel >= 0 ? (kdel < mdel ? kdel : mdel - framesize) : 0;
    //     for (n=0; n<offset; n++) {
    //     Complex   *fsigOut_ = (Complex *) fsigOut_.frame +NB*n;
    //     for (i = 0; i < NB; i++) fsigOut_[i].re = fsigOut_[i].im = FL(0.0);
    //     }
    //     for (n=offset; n<nsmps; n++) {
    //     Complex   *fsigIn = (Complex *) fsigIn.frame.auxp +NB*n;
    //     Complex   *fsigOut_ = (Complex *) fsigOut_.frame.auxp +NB*n;
    //     Complex   *delay = (Complex *) delframes.auxp +NB*n;

    //     for (i = 0; i < NB; i++) {
    //         delay[countr + i] = fsigIn[i];
    //         if (kdel) {
    //         if ((first = countr - kdel) < 0)
    //             first += mdel;

    //         for (j = first; j != countr; j = (j + framesize) % mdel) {
    //             amp += delay[j + i].re;
    //             freq += delay[j + i].im;
    //         }

    //         fsigOut_[i].re = (MYFLT) (amp / delayframes);
    //         fsigOut_[i].im = (MYFLT) (freq / delayframes);
    //         amp = freq = FL(0.0);
    //         }
    //         else {
    //         fsigOut_[i] = fsigIn[i];
    //         }
    //     }
    //     }
    //     countr += (N + 2);
    //     count = countr < mdel ? countr : 0;
    //     return OK;
    // }
    if(lastframe < fsigIn.framecount)
    {
        kdel = kdel >= 0 ? (kdel < mdel ? kdel : mdel - framesize) : 0;

        for(i = 0; i < N + 2; i += 2)
        {
            delay[countr + i]     = fin[i];
            delay[countr + i + 1] = fin[i + 1];

            if(kdel)
            {
                if((first = countr - kdel) < 0)
                    first += mdel;

                for(j = first; j != countr; j = (j + framesize) % mdel)
                {
                    amp += delay[j + i];
                    freq += delay[j + i + 1];
                }

                fout[i]     = (amp / delayframes);
                fout[i + 1] = (freq / delayframes);
                amp = freq = 0.;
            }
            else
            {
                fout[i]     = fin[i];
                fout[i + 1] = fin[i + 1];
            }
        }

        fsigOut_.framecount = lastframe = fsigIn.framecount;
        countr += (N + 2);
        count = countr < mdel ? countr : 0;
    }

    // return OK;
    // err1:
    // return csound->PerfError(csound, &(h),
    //                         Str("pvsblur: not initialised"));
    return fsigOut_;
}