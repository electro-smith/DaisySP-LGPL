#pragma once
#ifndef DSY_SPECTRALSCALE_H
#define DSY_SPECTRALSCALE_H

#include <cstddef>
#include "spectral.h"
#include "shy_fft.h"

namespace daicsp
{
/** Formant preservation methods for `SpectralScale`
 *  \param NONE - formants are not preserved
 *  \param LIFTERED - formants are preserved using a liftered cepstrum method
 *  \param ENVELOPE - formants are preserved using a true envelope method
 */
enum FORMANT
{
    NONE     = 0,
    LIFTERED = 1,
    ENVELOPE = 2,
};

/** SpectralScale
 * 
 *  Scales the pitch of a frequency-domain signal. 
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-06-14
 * 
 *  Ported from Csound pvscale
 */
template <size_t FFT_SIZE    = 2048,
          size_t OVERLAP     = FFT_SIZE / 4,
          size_t WINDOW_SIZE = FFT_SIZE>
class SpectralScale
{
  public:
    SpectralScale() {}
    ~SpectralScale() {}

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

    /** Initializes the SpectralFreezeFifo module.
             *  \param fsig_in - Initialized frequency-domain signal from the intended source.
             *  \param scale - Frequency multiplier.
             *  \param sample_rate - The program sample rate.
             *  \param formants - Formant preservation method. Defaults to none.
             *  \param gain - Amplitude scaling. Defaults to 1.
             *  \param coefficients - The number of cepstrum coefficients used in formant preservation. Defaults to 80.
             */
    void Init(SpectralBuffer<FFT_SIZE>& fsig_in,
              float                     scale,
              int                       sample_rate,
              FORMANT                   formants     = FORMANT::NONE,
              float                     gain         = 1.0f,
              int                       coefficients = 80);

    /** Processes an incoming `SpectralBuffer`
             *  \param fsig_in - An analyzed `SpectralBuffer`.
             *  \returns - A processed `SpectralBuffer`.
             */
    SpectralBuffer<FFT_SIZE>& Process(SpectralBuffer<FFT_SIZE>& fsig_in);

    void SetScale(float scale) { kscal_ = scale; }
    void SetFormants(FORMANT formant) { keepform_ = formant; }
    void SetGain(float gain) { gain_ = gain; }
    void SetCoefficients(float coefficients) { coefs_ = coefficients; }

    SpectralBuffer<FFT_SIZE>& GetFsig() { return fsig_out_; }

    STATUS GetStatus() { return status_; }

    void HaltOnError()
    {
        if(status_ != STATUS::OK)
            while(1) {}
    }

  private:
    int    sample_rate_;
    STATUS status_;

    SpectralBuffer<FFT_SIZE> fsig_out_;
    ShyFFT<float, FFT_SIZE>  fft_;

    float    kscal_;
    int      keepform_;
    float    gain_;
    int      coefs_;
    float    fenv_[FFT_SIZE + 2];
    float    ceps_[FFT_SIZE + 2];
    float    ceps_out_[FFT_SIZE + 2];
    float*   ftmp_;
    uint32_t lastframe;
};

#include "spectralscale.tcc"

} // namespace daicsp

#endif // DSY_SPECTRALSCALE_H
