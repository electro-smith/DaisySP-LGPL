#pragma once
#ifndef DSY_SPECTRALBLURFIFO_H
#define DSY_SPECTRALBLURFIFO_H

#include <cstddef>
#include "spectral.h"

namespace daicsp
{
/** SpectralBlurFifo
 * 
 *  Blurs a `SpectralBuffer` by smoothing frequency and time functions.
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-06-09
 *  Ported from Csound pvsblur
 */

template <size_t FFT_SIZE    = 2048,
          size_t OVERLAP     = 512,
          size_t WINDOW_SIZE = 2048>
class SpectralBlurFifo
{
  public:
    SpectralBlurFifo() {}
    ~SpectralBlurFifo() {}

    enum STATUS
    {
        OK = 0,
        E_SLIDING_NOT_IMPLEMENTED,
        E_FSIG_EQUAL,
        E_BUFFER_TOO_SMALL,
    };


    // void Init(SpectralBuffer& fsigIn, float delay, float maxDelay, int sampleRate, int block);

    /** Initializes the SpectralBlurFifo module.
         *  \param p - description
         */
    void Init(SpectralBuffer<FFT_SIZE + 2>& fsigIn,
              float                         delay,
              float*                        delayBuffer,
              size_t                        delaySize,
              int                           sampleRate,
              int                           block);

    /** Processes a single sample and returns it.
         *  \param in - input sample
         */
    SpectralBuffer<FFT_SIZE + 2>& Process(SpectralBuffer<FFT_SIZE + 2>& fsigIn);

    SpectralBuffer<FFT_SIZE + 2>& GetFsig() { return fsigOut_; }

    STATUS GetStatus() { return status_; }

    /** Enters an infinite loop if the module encountered an error.
         *  Useful for error checking.
         */
    void HaltOnError()
    {
        if(status_ != STATUS::OK)
            while(1)
                ;
    }

  private:
    SpectralBuffer<FFT_SIZE + 2> fsigOut_;
    STATUS                       status_;
    int                          sr_;

    float        kdel_;
    float        maxdel;
    float*       delframes;
    float        frpsec;
    int          count;
    unsigned int lastframe;
};

#include "spectralblurfifoimpl.h"

} // namespace daicsp

#endif // DSY_SPECTRALBLURFIFO_H
