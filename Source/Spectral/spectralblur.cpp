#include <cstddef>
#include "spectral.h"

#include "spectralblur.h"

using namespace daicsp;

void SpectralBlur::Init(SpectralBuffer& fsig_in,
                        float           delay,
                        float*          delay_buffer,
                        size_t          delay_size,
                        int             sample_rate)
{
    status_      = STATUS::OK;
    sample_rate_ = sample_rate;
    kdel_        = delay;
    delframes    = delay_buffer;

    float* delay_ptr;
    int    N    = fsig_in.N, i, j;
    int    olap = fsig_in.overlap;
    int    delayframes, framesize = N + 2;

    // Calculating the max delay time from the given buffer
    maxdel = (float)delay_size / ((N + 2) * ((float)sample_rate / olap));

    if(maxdel < delay)
    {
        status_ = STATUS::E_BUFFER_TOO_SMALL;
        return;
    }

    //   if (fsigIn == fsig_out_)
    //     csound->Warning(csound, Str("Unsafe to have same fsig as in and out"));
    if(&fsig_in == &fsig_out_)
    {
        status_ = STATUS::E_FSIG_EQUAL;
        return;
    }

    if(fsig_in.sliding)
    {
        status_ = STATUS::E_SLIDING_NOT_IMPLEMENTED;
        return;
        // csound->InitError(csound, Str("pvsblur does not work sliding yet"));
        // delayframes = (int) (FL(0.5) + *maxdel * CS_ESR);
        // if (fsig_out_.frame.auxp == NULL ||
        //     fsig_out_.frame.size < sizeof(MYFLT) * CS_KSMPS * (N + 2))
        //   csound->AuxAlloc(csound, (N + 2) * sizeof(MYFLT) * CS_KSMPS,
        //                    &fsig_out_.frame);

        // if (delframes.auxp == NULL ||
        //     delframes.size < (N + 2) * sizeof(MYFLT) * CS_KSMPS * delayframes)
        //   csound->AuxAlloc(csound,
        //                    (N + 2) * sizeof(MYFLT) * CS_KSMPS * delayframes,
        //                    &delframes);
    }
    else
    {
        frpsec      = (float)sample_rate_ / olap;
        delayframes = (int)(maxdel * frpsec);

        //   NOTE -- No need for allocation with this scheme
        //   if (fsig_out_.frame.auxp == NULL ||
        //       fsig_out_.frame.size < sizeof(float) * (N + 2))
        //     csound->AuxAlloc(csound, (N + 2) * sizeof(float), &fsig_out_.frame);

        //   if (delframes.auxp == NULL ||
        //       delframes.size < (N + 2) * sizeof(float) * CS_KSMPS * delayframes)
        //     csound->AuxAlloc(csound, (N + 2) * sizeof(float) * delayframes,
        //                      &delframes);
    }
    delay_ptr = delframes;

    for(j = 0; j < framesize * delayframes; j += framesize)
        for(i = 0; i < N + 2; i += 2)
        {
            delay_ptr[i + j]     = 0.0f;
            delay_ptr[i + j + 1] = i * (float)sample_rate_ / N;
        }

    fsig_out_.N          = N;
    fsig_out_.overlap    = olap;
    fsig_out_.winsize    = fsig_in.winsize;
    fsig_out_.wintype    = fsig_in.wintype;
    fsig_out_.format     = fsig_in.format;
    fsig_out_.framecount = 1;
    lastframe            = 0;
    count                = 0;
    fsig_out_.sliding    = fsig_in.sliding;
    fsig_out_.NB         = fsig_in.NB;
    fsig_out_.ready      = false;
}

SpectralBuffer& SpectralBlur::Process(SpectralBuffer& fsig_in)
{
    if(fsig_in.ready)
    {
        int    j, i, N = fsig_out_.N, first, framesize = N + 2;
        int    countr = count;
        float  amp = 0.0, freq = 0.0;
        int    delayframes = (int)(kdel_ * frpsec);
        int    kdel        = delayframes * framesize;
        int    mdel        = (int)(maxdel * frpsec) * framesize;
        float* fin         = fsig_in.frame;
        float* fout        = fsig_out_.frame;
        float* delay       = delframes;

        // if (UNLIKELY(fsig_out_ == NULL || delay == NULL)) goto err1;

        // NOTE -- sliding doesn't work yet
        // if (fsigIn.sliding) {
        //     unsigned int offset = 0;
        //     unsigned int n, nsmps = block;
        //     int NB = fsigIn.NB;
        //     kdel = kdel >= 0 ? (kdel < mdel ? kdel : mdel - framesize) : 0;
        //     for (n=0; n<offset; n++) {
        //     Complex   *fsig_out_ = (Complex *) fsig_out_.frame +NB*n;
        //     for (i = 0; i < NB; i++) fsig_out_[i].re = fsig_out_[i].im = FL(0.0);
        //     }
        //     for (n=offset; n<nsmps; n++) {
        //     Complex   *fsigIn = (Complex *) fsigIn.frame.auxp +NB*n;
        //     Complex   *fsig_out_ = (Complex *) fsig_out_.frame.auxp +NB*n;
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

        //         fsig_out_[i].re = (MYFLT) (amp / delayframes);
        //         fsig_out_[i].im = (MYFLT) (freq / delayframes);
        //         amp = freq = FL(0.0);
        //         }
        //         else {
        //         fsig_out_[i] = fsigIn[i];
        //         }
        //     }
        //     }
        //     countr += (N + 2);
        //     count = countr < mdel ? countr : 0;
        //     return OK;
        // }
        if(lastframe < fsig_in.framecount)
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

            fsig_out_.framecount = lastframe = fsig_in.framecount;
            countr += (N + 2);
            count = countr < mdel ? countr : 0;
        }
    }

    // return OK;
    // err1:
    // return csound->PerfError(csound, &(h),
    //                         Str("pvsblur: not initialised"));
    fsig_out_.ready = fsig_in.ready;
    return fsig_out_;
}