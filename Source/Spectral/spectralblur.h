#pragma once
#ifndef DSY_SPECTRALBLUR_H
#define DSY_SPECTRALBLUR_H

#include <cstddef>

namespace daicsp
{
/** SpectralBlur
 * 
 *  Blurs a `SpectralBuffer` by smoothing frequency and time functions.
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-06-09
 *  Ported from Csound pvsblur
 */

class SpectralBlur
{
  public:
    SpectralBlur() {}
    ~SpectralBlur() {}

    enum STATUS
    {
        OK = 0,
        E_SLIDING_NOT_IMPLEMENTED,
        E_FSIG_EQUAL,
        E_BUFFER_TOO_SMALL,
    };


    // void Init(SpectralBuffer& fsigIn, float delay, float maxDelay, int sampleRate, int block);

    /** Initializes the SpectralBlur module.
         *  \param p - description
         */
    void Init(SpectralBuffer<FFT::MAX_FLOATS>& fsigIn,
              float                            delay,
              float*                           delayBuffer,
              size_t                           delaySize,
              int                              sampleRate,
              int                              block);

    /** Processes a single sample and returns it.
         *  \param in - input sample
         */
    SpectralBuffer<FFT::MAX_FLOATS>&
    Process(SpectralBuffer<FFT::MAX_FLOATS>& fsigIn, int block);

    SpectralBuffer<FFT::MAX_FLOATS>& GetFsig() { return fsigOut_; }

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
    SpectralBuffer<FFT::MAX_FLOATS> fsigOut_;
    STATUS                          status_;
    int                             sr_;

    float        kdel_;
    float        maxdel;
    float*       delframes;
    float        frpsec;
    int          count;
    unsigned int lastframe;
};

} // namespace daicsp

#endif // DSY_SPECTRALBLUR_H
