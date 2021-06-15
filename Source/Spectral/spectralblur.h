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
class SpectralBlur
{
  public:
    SpectralBlur() {}
    ~SpectralBlur() {}


    /** Indicates the current status of the module. 
         *  Warnings are indicated by a leading W, and are silently corrected. 
         *  Errors are indicated by a leading E and cause an immediate exit.
         * 
         *  \param OK - No errors have been reported.
         *  \param E_FSIG_EQUAL - The input and output `SpectralBuffer` cannot be the same.
         *  \param E_BUFFER_TOO_SMALL - The given buffer is too small for the specified delay time.
         *  \param E_SLIDING_NOT_IMPLEMENTED - Sliding is currently not implemented, so the overlap size must be greater than the audio block size.
         */
    enum STATUS
    {
        OK = 0,
        E_FSIG_EQUAL,
        E_BUFFER_TOO_SMALL,
        E_SLIDING_NOT_IMPLEMENTED,
    };

    // void Init(SpectralBuffer& fsigIn, float delay, float maxDelay, int sampleRate, int block);

    /** Initializes the SpectralBlurFifo module.
         *  \param fsig_in - Initialized frequency-domain signal from the intended source.
         *  \param delay - The delay time in seconds. The supplied buffer must be large enough
         *  to hold this amount of delay. Each second of delay corresponds to 
         *  `(FFT_SIZE + 2) * (sample_rate / OVERLAP)` number of floats.
         *  \param delay_buffer - Pointer to a buffer located in the Daisy's SDRAM BSS section.
         *  This can be declared as `float DSY_SDRAM_BSS myBuff[MY_BUFF_SIZE];`
         *  \param delay_size - The size of the given delay buffer.
         *  \param sample_rate - The program sample rate. 
         */
    void Init(SpectralBuffer& fsig_in,
              float           delay,
              float*          delay_buffer,
              size_t          delay_size,
              int             sample_rate);

    /** Processes an incoming `SpectralBuffer`
         *  \param fsig_in - An analyzed `SpectralBuffer`.
         *  \returns - A processed `SpectralBuffer`.
         */
    SpectralBuffer& Process(SpectralBuffer& fsig_in);

    void SetDelay(float delay) { kdel_ = delay; }

    SpectralBuffer& GetFsig() { return fsig_out_; }

    STATUS GetStatus() { return status_; }

    /** Enters an infinite loop if the module encountered an error.
         *  Useful for error checking.
         */
    void HaltOnError()
    {
        if(status_ != STATUS::OK)
            while(1) {}
    }

  private:
    SpectralBuffer fsig_out_;
    STATUS         status_;
    int            sample_rate_;

    float        kdel_;
    float        maxdel;
    float*       delframes;
    float        frpsec;
    int          count;
    unsigned int lastframe;
};

} // namespace daicsp

#endif // DSY_SPECTRALBLURFIFO_H
