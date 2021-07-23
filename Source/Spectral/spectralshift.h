#pragma once
#ifndef DSY_SPECTRALSHIFT_H
#define DSY_SPECTRALSHIFT_H

#include <cstddef>
#include <cstdint>
#include "spectral.h"

namespace daicsp
{

/** SpectralShift
 * 
 *  Shift the pitch of a frequency-domain signal.
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-07-23
 * 
 *  Ported from Csound pvshift
 */
class SpectralShift
{
  public:
    SpectralShift() {}
    ~SpectralShift() {}

    /** Indicates the current status of the module. 
             *  Warnings are indicated by a leading W, and are silently corrected. 
             *  Errors are indicated by a leading E and cause an immediate exit.
             * 
             *  \param OK - No errors have been reported.
             *  \param E_FSIG_EQUAL - The input and output `SpectralBuffer` cannot be the same.
             *  \param W_KEEPFORM_NOT_IMPLEMENTED - Any keepform mode other than NONE is not yet supported.
             *  \param E_SLIDING_NOT_IMPLEMENTED - Sliding is currently not implemented, so the overlap size must be greater than the audio block size.
             */
    enum STATUS
    {
        OK = 0,
        E_FSIG_EQUAL,
        W_KEEPFORM_NOT_IMPLEMENTED,
        E_SLIDING_NOT_IMPLEMENTED,
    };

    /** Initializes the SpectralShift module.
             *  \param fsig_in Initialized frequency-domain signal from the intended source.
             *  \param shift Shift amount in Hz (can be positive or negative).
             *  \param lowest Lowest frequency to shift.
             *  \param sample_rate The program sample rate.
             *  \param keepform Formant preservation mode. Defaults to none. (other modes not yet supported)
             *  \param gain Aplitude scaling (defaults to 1).
             *  \param coefficients The number of cepstrum coefficients used in formant preservation (defaults to 80).
             */
    void Init(SpectralBuffer& fsig_in,
              float           shift,
              float           lowest,
              int             sample_rate,
              FORMANT         keepform = FORMANT::NONE,
              float           gain = 1,
              int             coefficients = 80);

    /** Processes an incoming `SpectralBuffer`
             *  \param fsig_in - An analyzed `SpectralBuffer`.
             *  \returns - A processed `SpectralBuffer`.
             */
    SpectralBuffer& ParallelProcess(SpectralBuffer& fsig_in);

    void SetShift(float shift) { shift_ = shift; }
    void SetLowest(float lowest) { lowest_ = lowest; }
    void SetKeepform(FORMANT keepform) { keepform_ = keepform; }
    void SetGain(float gain) { gain_ = gain; }
    void SetCoefficients(int coefficients) { coefficients_ = coefficients; }

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

    SpectralBuffer fsig_out_;

    float    fenv_[kFFTMaxFloats];
    // float    ceps_[kFFTMaxFloats]; // NOTE -- this can be removed if we don't expect to use keepform modes
    float    ftmp_[kFFTMaxFloats];

    float    shift_; 
    float    lowest_;
    FORMANT  keepform_;
    float    gain_; 
    int      coefficients_;
    uint32_t lastframe_;
};

} // namespace daicsp

#endif // DSY_SPECTRALSHIFT_H
