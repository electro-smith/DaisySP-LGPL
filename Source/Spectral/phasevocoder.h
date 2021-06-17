#pragma once
#ifndef DSY_PHASEVOCODERFIFO_H
#define DSY_PHASEVOCODERFIFO_H

#include <cstring>
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
class PhaseVocoder
{
  public:
    PhaseVocoder() {}
    ~PhaseVocoder() {}

    /** Indicates the current status of the module. 
         *  Warnings are indicated by a leading W, and are silently corrected. 
         *  Errors are indicated by a leading E and cause an immediate exit.
         * 
         *  \param OK - No errors have been reported.
         *  \param E_BLOCK_TOO_BIG - The audio block in combination with the fft size is too large for the input buffer. 
         *  \param W_BUFFER_UNDERFLOW - The input buffer was filled before the previous buffer was fully processed.
         *  \param W_BUFFER_MISMATCH - The incoming signal is ready before any output is requested.
         *  \param W_INVALID_STATE - Against all odds, you've put the state_ property into an invalid state. 
         *  \param E_SLIDING_NOT_IMPLEMENTED - Sliding is currently not implemented, so the overlap size must be greater than the audio block size.
         */
    enum STATUS
    {
        OK = 0,
        E_BLOCK_TOO_BIG,
        W_BUFFER_UNDERFLOW,
        W_BUFFER_MISMATCH,
        W_INVALID_STATE,
        E_SLIDING_NOT_IMPLEMENTED,
    };

    /** Initializes the PhaseVocoderFifo module.
         *  \param fsig_in - Initialized frequency-domain signal from the intended source.
         *  \param sample_rate - The program sample rate.
         */
    void Init(SpectralBuffer& fsig_in, size_t sample_rate, size_t audio_block);

    /** Processes an incoming fsig in parallel to the audio callback.
         *  \param fsig_in - A `SpectralBuffer` from the source used to initialize this module.
         */
    void ParallelProcess(SpectralBuffer& fsig_in); // pvsynth

    /** Retrieves a single sample from the FIFO, and
     *  queues the bulk processing when appropriate.
     */
    float Process();

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

    /** Check if the vocoder is ready for processing
     * 
     */
    int GetQueuedFrameCount() { return output_count_; }

  private:
    /** Corresponds to pvsynth's pvssynthset -- Phase Vocoder Synthesis _sliding_ set.
         *  This is not currently implemented, but can be useful for small overlap sizes.
         */
    void ProcessSliding(SpectralBuffer& fsig_in,
                        size_t          size); // pvssynth

    float Tick(SpectralBuffer& fsig_in); // analyze_tick

    void GenerateFrame(SpectralBuffer& fsig_in); // process_frame

    int    buflen_;
    float  RoverTwoPi_, TwoPioverR_, Fexact_;
    float* nextOut_;
    int    nO_, Ii_, IOi_;
    size_t outptr_;

    float  output_[kFFTMaxFrames];
    float  overlapbuf_[kFFTMaxOverlapBuff];
    size_t overlap_;
    size_t num_overlaps_;
    size_t output_count_;
    float* output_segment_;
    float* process_segment_;

    float synbuf_[kFFTMaxFloats];
    float synbufOut_[kFFTMaxFloats];
    float analwinbuf_[kFFTMaxWindow];
    float synwinbuf_[kFFTMaxWindow];
    float oldOutPhase_[kFFTMaxBins];

    ShyFFT<float, kFFTMaxSize> fft_;

    float sample_rate_;

    /* check these against fsig vals */
    // int overlap,winsize,fftsize,wintype,format;
    int    bin_index_; /* for phase normalization across frames */
    STATUS status_;
};

} // namespace daicsp

#endif // DSY_PHASEVOCODERFIFO_H
