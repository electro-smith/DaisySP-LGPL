#pragma once
#ifndef DSY_SPECTRALSMOOTH_H
#define DSY_SPECTRALSMOOTH_H

#include <cstddef>
#include <cstdint>
#include "spectral.h"

namespace daicsp
{

/** SpectralSmooth
 * 
 *  Smooths the input frequency spectrum over time. 
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-07-23
 * 
 *  Ported from Csound pvsmooth
 */
class SpectralSmooth
{
  public:
    SpectralSmooth() {}
    ~SpectralSmooth() {}

    /** Indicates the current status of the module. 
             *  Warnings are indicated by a leading W, and are silently corrected. 
             *  Errors are indicated by a leading E and cause an immediate exit.
             * 
             *  \param OK - No errors have been reported.
             *  \param E_FSIG_EQUAL - The input and output `SpectralBuffer` cannot be the same.
             *  \param E_SLIDING_NOT_IMPLEMENTED - Sliding is currently not implemented, so the overlap size must be greater than the audio block size.
             */
    enum STATUS
    {
        OK = 0,
        E_FSIG_EQUAL,
        E_SLIDING_NOT_IMPLEMENTED,
    };

    /** Initializes the SpectralSmooth module.
             *  \param fsig_in Initialized frequency-domain signal from the intended source.
             *  \param scale Amount of cutoff frequency for amplitude function filtering (0 to 1).
             *  \param kfcf Amount of cutoff frequency for frequency function filtering (0 to 1).
             *  \param sample_rate The program sample rate.
             */
    void Init(SpectralBuffer& fsig_in,
              float           kacf,
              float           kfcf,
              int             sample_rate);

    /** Processes an incoming `SpectralBuffer`
             *  \param fsig_in - An analyzed `SpectralBuffer`.
             *  \returns - A processed `SpectralBuffer`.
             */
    SpectralBuffer& ParallelProcess(SpectralBuffer& fsig_in);

    /** Sets the amount of cutoff frequency for amplitude function filtering
     *  \param kacf amount between 0 and 1 (clamped)
     */
    void SetKacf(float kacf) { kfra_ = kacf; }

    /** Sets the amount of cutoff frequency for frequency function filtering
     *  \param kfcf amount between 0 and 1 (clamped)
     */
    void SetKfcf(float kfcf) { kfrf_ = kfcf; }

    SpectralBuffer& GetFsig() { return fsig_out_; }

    STATUS GetStatus() { return status_; }

    void HaltOnError()
    {
        if(status_ != STATUS::OK)
            while(1) {}
    }

  private:
    int    sample_rate_;
    STATUS status_;

    SpectralBuffer             fsig_out_;

    // NOTE -- these two names do not reflect their official documentation
    float    kfra_; // -> kacf
    float    kfrf_; // -> kfcf
    float    del_[kFFTMaxFloats];
    uint32_t lastframe_;
};

} // namespace daicsp

#endif // DSY_SPECTRALSMOOTH_H
