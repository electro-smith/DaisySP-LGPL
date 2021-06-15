#include <cstring>

#include "spectral.h"
#include "spectralfreeze.h"

using namespace daicsp;

void SpectralFreeze::Init(SpectralBuffer& fsig_in,
                          float           freeze_amp,
                          float           freeze_freq,
                          int             sample_rate)
{
    sample_rate_ = sample_rate;
    status_      = STATUS::OK;
    int N        = fsig_in.N;
    kfra_        = freeze_amp;
    kfrf_        = freeze_freq;

    fsig_out_.N          = N;
    fsig_out_.overlap    = fsig_in.overlap;
    fsig_out_.winsize    = fsig_in.winsize;
    fsig_out_.wintype    = fsig_in.wintype;
    fsig_out_.format     = fsig_in.format;
    fsig_out_.framecount = 1;
    fsig_out_.ready      = false;
    lastframe_           = 0;

    fsig_out_.NB      = (N / 2) + 1;
    fsig_out_.sliding = fsig_in.sliding;
    if(fsig_in.sliding)
    {
        status_ = STATUS::E_SLIDING_NOT_IMPLEMENTED;
        return;

        // uint32_t nsmps = CS_KSMPS;
        // if (fsigOut_.frame.auxp == NULL ||
        //     fsigOut_.frame.size < sizeof(MYFLT) * (N + 2) * nsmps)
        //   csound->AuxAlloc(csound, (N + 2) * sizeof(MYFLT) * nsmps,
        //                    &fsigOut_.frame);
        // if (freez.auxp == NULL ||
        //     freez.size < sizeof(MYFLT) * (N + 2) * nsmps)
        //   csound->AuxAlloc(csound, (N + 2) * sizeof(MYFLT) * nsmps, &freez);
    }
    else
    {
        // if (fsigOut_.frame.auxp == NULL ||
        //     fsigOut_.frame.size < sizeof(float) * (N + 2))
        //   csound->AuxAlloc(csound, (N + 2) * sizeof(float), &fsigOut_.frame);
        // if (freez.auxp == NULL || freez.size < sizeof(float) * (N + 2))
        //   csound->AuxAlloc(csound, (N + 2) * sizeof(float), &freez);

        // if (UNLIKELY(!((fsigOut_.format == PVS_AMP_FREQ) ||
        //                (fsigOut_.format == PVS_AMP_PHASE))))
        //   return csound->InitError(csound, Str("pvsfreeze: signal format "
        //                                        "must be amp-phase or amp-freq."));
    }
    //   return OK;
}

SpectralBuffer& SpectralFreeze::Process(SpectralBuffer& fsig_in)
{
    if(fsig_in.ready)
    {
        int    i;
        int    framesize;
        float  freeza, freezf;
        float *fout, *fin, *freez;
        // NOTE -- sliding not yet implemented
        // if (fsigIn.sliding)
        //   return pvssfreezeprocess(csound, p);

        // Preventing wayward changes if interrupt occurs between assignments
        {
            // TODO -- add irq blocker to prevent parameter mismatch
            // daisy::ScopedIrqBlocker block;
            freeza = kfra_;
            freezf = kfrf_;
        }


        fout  = fsig_out_.frame;
        fin   = fsig_in.frame;
        freez = freez_;
        int N = fsig_in.N;

        framesize = N + 2;

        if(lastframe_ < fsig_in.framecount)
        {
            memset(fout, 0, sizeof(float) * (N + 2));
            for(i = 0; i < framesize; i += 2)
            {
                if(freeza < 1)
                    freez[i] = fin[i];
                if(freezf < 1)
                    freez[i + 1] = fin[i + 1];
                fout[i]     = freez[i];
                fout[i + 1] = freez[i + 1];
            }
            fsig_out_.framecount = lastframe_ = fsig_in.framecount;
        }
    }

    fsig_out_.ready = fsig_in.ready;
    return fsig_out_;
}