#pragma once
#ifndef DSY_PHASEVOCODERFIFO_H
#define DSY_PHASEVOCODERFIFO_H

#include <math.h>
#include <string.h>
#include "spectral.h"
#include "shy_fft.h"

namespace daicsp
{
/** PhaseVocoder
 * 
 *  Converts frequency-domain signals back to time-domain using a FIFO interface.
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-06-03
 *  Ported from Csound pvsynth.
 */
template <size_t FFT_SIZE    = 2048,
          size_t OVERLAP     = FFT_SIZE / 4,
          size_t WINDOW_SIZE = FFT_SIZE>
class PhaseVocoderFifo
{
  public:
    PhaseVocoderFifo() {}
    ~PhaseVocoderFifo() {}

    /** Indicates the current status of the module. 
         *  Warnings are indicated by a leading W, and are silently corrected. 
         *  Errors are indicated by a leading E and cause an immediate exit.
         * 
         *  \param OK - No errors have been reported.
         *  \param W_BUFFER_UNDERFLOW - The input buffer was filled before the previous buffer was fully processed.
         *  \param W_INVALID_STATE - Against all odds, you've put the state_ property into an invalid state. 
         *  \param E_SLIDING_NOT_IMPLEMENTED - Sliding is currently not implemented, so the overlap size must be greater than the audio block size.
         */
    enum STATUS
    {
        OK = 0,
        E_SLIDING_NOT_IMPLEMENTED,
        W_BUFFER_UNDERFLOW,
        W_INVALID_STATE,
    };

    enum STATE
    {
        INIT = 0,
        IDLE,
        PROCESSING,
    };

    /** Initializes the PhaseVocoderFifo module.
         *  \param fsig_in - Initialized frequency-domain signal from the intended source.
         *  \param sample_rate - The program sample rate.
         */
    void Init(SpectralBuffer<FFT_SIZE>& fsig_in, size_t sample_rate);

    /** Processes an incoming fsig.
         *  \param fsig_in - A `SpectralBuffer` from the source used to initialize this module.
         */
    void Process(SpectralBuffer<FFT_SIZE>& fsig_in); // pvsynth

    /** Retrieves a single sample from the FIFO, and
     *  queues the bulk processing when appropriate.
     */
    float Sample();

    /** Retrieves the current status. Useful for error checking.
         */
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
    /** Corresponds to pvsynth's pvssynthset -- Phase Vocoder Synthesis _sliding_ set.
         *  This is not currently implemented, but can be useful for small overlap sizes.
         */
    void ProcessSliding(SpectralBuffer<FFT_SIZE>& fsig_in,
                        size_t                    size); // pvssynth

    float Tick(SpectralBuffer<FFT_SIZE>& fsig_in); // analyze_tick

    void GenerateFrame(SpectralBuffer<FFT_SIZE>& fsig_in); // process_frame

    int    buflen_;
    float  RoverTwoPi_, TwoPioverR_, Fexact_;
    float* nextOut_;
    int    nO_, Ii_, IOi_;
    size_t outptr_;

    float  output_[FFT_SIZE * 4];
    float  overlapbuf_[OVERLAP * 2];
    size_t half_overlap_;
    float* output_segment_;
    float* process_segment_;

    float synbuf_[FFT_SIZE + 2];
    float synbufOut_[FFT_SIZE + 2];
    float analwinbuf_[FFT_SIZE + 1];
    float synwinbuf_[FFT_SIZE + 1];
    float oldOutPhase_[FFT_SIZE / 2 + 1];

    ShyFFT<float, FFT_SIZE> fft_;

    float sample_rate_;

    /* check these against fsig vals */
    // int overlap,winsize,fftsize,wintype,format;
    int    bin_index_; /* for phase normalization across frames */
    STATUS status_;
    STATE  state_;
};

#include "phasevocoderfifo.tcc"

} // namespace daicsp

#endif // DSY_PHASEVOCODERFIFO_H
