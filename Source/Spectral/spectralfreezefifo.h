#pragma once
#ifndef DSY_SPECTRALFREEZEFIFO_H
#define DSY_SPECTRALFREEZEFIFO_H

#include <cstddef>

namespace daicsp
{
/** SpectralFreezeFifo
 *  
 *  Author: Gabriel Ball
 *  Date: 2021-06-11
 *  
 *  Ported from Csound pvsfreeze
 */
template <size_t FFT_SIZE    = 2048,
          size_t OVERLAP     = 512,
          size_t WINDOW_SIZE = 2048>
class SpectralFreezeFifo
{
  public:
    SpectralFreezeFifo() {}
    ~SpectralFreezeFifo() {}

    /** Indicates the current status of the module. 
         *  Warnings are indicated by a leading W, and are silently corrected. 
         *  Errors are indicated by a leading E and cause an immediate exit.
         * 
         *  \param OK - No errors have been reported.
         *  \param E_SLIDING_NOT_IMPLEMENTED - Sliding is currently not implemented, so the overlap size must be greater than the audio block size.
         */
    enum STATUS
    {
        OK = 0,
        E_SLIDING_NOT_IMPLEMENTED,
    };

    /** Initializes the SpectralFreezeFifo module.
         *  \param fsig_in - Initialized frequency-domain signal from the intended source.
         *  \param freeze_amp - Freezing switch for amplitudes. On if above or equal to 1, and off otherwise.
         *  \param freeze_freq - Freezing switch for frequencies. On if above or equal to 1, and off otherwise.
         *  \param sample_rate - The program sample rate.
         */
    void Init(SpectralBuffer<FFT_SIZE>& fsig_in,
              float                     freeze_amp,
              float                     freeze_freq,
              int                       sample_rate);

    /** Processes an incoming `SpectralBuffer`
         *  \param fsig_in - An analyzed `SpectralBuffer`.
         *  \returns - A processed `SpectralBuffer`.
         */
    SpectralBuffer<FFT_SIZE>& Process(SpectralBuffer<FFT_SIZE>& fsig_in);

    void SetAmplitude(float amp) { kfra_ = amp; }

    void SetFrequency(float freq) { kfrf_ = freq; }

    STATUS GetStatus() { return status_; }

    void HaltOnError()
    {
        if(status_ != STATUS::OK)
            while(1)
                ;
    }


  private:
    int    sample_rate_;
    STATUS status_;

    SpectralBuffer<FFT_SIZE> fsig_out_;

    float    kfra_, kfrf_;
    float    freez_[FFT_SIZE + 2];
    uint32_t lastframe_;
};

#include "spectralfreezefifo.tcc"

} // namespace daicsp

#endif // DSY_SPECTRALFREEZEFIFO_H
