#include <cmath>
#include <cstring>
#include "spectralsmooth.h"
#include "dsy_lgpl_dsp.h"


using namespace daicsp;

void SpectralSmooth::Init(SpectralBuffer& fsig_in,
              float           kacf,
              float           kfcf,
              int             sample_rate)
{
  status_      = STATUS::OK;
  sample_rate_ = sample_rate;
  kfra_ = kacf;
  kfrf_ = kfcf;

  int    N = fsig_in.N;

  if(fsig_in == fsig_out_)
  {
      status_ = STATUS::E_FSIG_EQUAL;
      return;
  }

  fsig_out_.NB = (N/2)+1;
  fsig_out_.sliding = fsig_in.sliding;
  if (fsig_in.sliding) {
    // if (p->fout->frame.auxp == NULL ||
    //     p->fout->frame.size < sizeof(MYFLT) * CS_KSMPS * (N + 2))
    //   csound->AuxAlloc(csound, (N + 2) * sizeof(MYFLT) * CS_KSMPS,
    //                    &p->fout->frame);
    // if (p->del.auxp == NULL ||
    //     p->del.size < sizeof(MYFLT) * CS_KSMPS * (N + 2))
    //   csound->AuxAlloc(csound, (N + 2) * sizeof(MYFLT) * CS_KSMPS,
    //                    &p->del);

    status_ = STATUS::E_SLIDING_NOT_IMPLEMENTED;
    return;
  }
  // else
  // {
  //   if (p->fout->frame.auxp == NULL ||
  //       p->fout->frame.size < sizeof(float) * (N + 2))
  //     csound->AuxAlloc(csound, (N + 2) * sizeof(float), &p->fout->frame);
  // NOTE -- the size of the delay buffer is just kMaxFloats
  //   if (p->del.auxp == NULL || p->del.size < sizeof(float) * (N + 2))
  //     csound->AuxAlloc(csound, (N + 2) * sizeof(float), &p->del);
  // }
  memset(del_, 0, (N + 2) * sizeof(float));
  fsig_out_.N = N;
  fsig_out_.overlap = fsig_in.overlap;
  fsig_out_.winsize = fsig_in.winsize;
  fsig_out_.wintype = fsig_in.wintype;
  fsig_out_.format = fsig_in.format;
  fsig_out_.framecount = 1;
  fsig_out_.ready = false;
  lastframe_ = 0;
  // NOTE -- no need to worry about this, since csound hasn't implemented anything else anyway
  // if (UNLIKELY(!((p->fout->format == PVS_AMP_FREQ) ||
  //                (p->fout->format == PVS_AMP_PHASE))))
  //   return csound->InitError(csound, Str("pvsmooth: signal format "
  //                                        "must be amp-phase or amp-freq."));
  // return OK;
}

SpectralBuffer& SpectralSmooth::ParallelProcess(SpectralBuffer& fsig_in)
{
  if (fsig_in.ready)
  {
    //  IGN(csound); -- what is this, an optional assertion or something?
    int32_t i;
    int     framesize;
    float  ffa, ffr;

    ffa = kfra_;
    ffr = kfrf_;

    framesize = fsig_in.N + 2;

    // NOTE -- sliding not yet implemented
    // if (p->fin->sliding) {
    //   CMPLX *fout, *fin, *del;
    //   double  costh1, costh2, coef1, coef2;
    //   uint32_t offset = p->h.insdshead->ksmps_offset;
    //   uint32_t n, nsmps = CS_KSMPS;
    //   int32_t NB = p->fin->NB;
    //   ffa = ffa < 0.0 ? 0.0 : (ffa > 1.0 ? 1.0 : ffa);
    //   ffr = ffr < 0.0 ? 0.0 : (ffr > 1.0 ? 1.0 : ffr);
    //   costh1 = 2.0 - cos(PI * ffa);
    //   costh2 = 2.0 - cos(PI * ffr);
    //   coef1 = sqrt(costh1 * costh1 - 1.0) - costh1;
    //   coef2 = sqrt(costh2 * costh2 - 1.0) - costh2;

    //   for (n=0; n<offset; n++)
    //     for (i=0; i<NB; i++) {
    //       fout = (CMPLX*) p->fout->frame.auxp +NB*n;
    //       del = (CMPLX*) p->del.auxp +NB*n;
    //       fout[i].re = fout[i].im = del[i].re = del[i].im = FL(0.0);
    //     }
    //   for (n=offset; n<nsmps; n++) {
    //     fout = (CMPLX*) p->fout->frame.auxp +NB*n;
    //     fin = (CMPLX*) p->fin->frame.auxp +NB*n;
    //     del = (CMPLX*) p->del.auxp +NB*n;
    //     if (IS_ASIG_ARG(p->kfra)) {
    //       ffa = (double)  p->kfra[n];
    //       ffa = ffa < 0.0 ? 0.0 : (ffa > 1.0 ? 1.0 : ffa);
    //       costh1 = 2.0 - cos(PI * ffa);
    //       coef1 = sqrt(costh1 * costh1 - 1.0) - costh1;
    //     }
    //     if (IS_ASIG_ARG(p->kfrf)) {
    //       ffr = (double)  p->kfrf[n];
    //       ffr = ffr < 0.0 ? 0.0 : (ffr > 1.0 ? 1.0 : ffr);
    //       costh2 = 2.0 - cos(PI * ffr);
    //       coef2 = sqrt(costh2 * costh2 - 1.0) - costh2;
    //     }
    //     for (i=0; i<NB; i++) {
    //       /* amp smoothing */
    //       fout[i].re = fin[i].re * (1.0 + coef1) - del[i].re * coef1;
    //       /* freq smoothing */
    //       fout[i].im = fin[i].im * (1.0 + coef2) - del[i].im * coef2;
    //       del[i] = fout[i];
    //     }
    //   }
    //   return OK;
    // }
    if (lastframe_ < fsig_in.framecount) {
      float   *fout, *fin, *del;
      float  costh1, costh2, coef1, coef2;
      fout = (float *) fsig_out_.frame;
      fin = (float *) fsig_in.frame;
      del = (float *) del_;

      ffa = ffa < 0.0f ? 0.0f : (ffa > 1.0f ? 1.0f : ffa);
      ffr = ffr < 0.0f ? 0.0f : (ffr > 1.0f ? 1.0f : ffr);
      costh1 = 2.0 - cos(PI_F * ffa);
      costh2 = 2.0 - cos(PI_F * ffr);
      coef1 = sqrt(costh1 * costh1 - 1.0f) - costh1;
      coef2 = sqrt(costh2 * costh2 - 1.0f) - costh2;

      for (i = 0; i < framesize; i += 2) {
        /* amp smoothing */
        fout[i] = (float) (fin[i] * (1.0 + coef1) - del[i] * coef1);
        /* freq smoothing */
        fout[i + 1] = (float) (fin[i + 1] * (1.0 + coef2) - del[i + 1] * coef2);
        del[i] = fout[i];
        del[i + 1] = fout[i + 1];
      }
      fsig_out_.framecount = lastframe_ = fsig_in.framecount;
    }
  }
  
  fsig_out_.ready = fsig_in.ready;
  return fsig_out_;
}
