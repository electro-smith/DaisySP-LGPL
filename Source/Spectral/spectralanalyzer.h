#pragma once
#ifndef DSY_SPECTRALANALYZERFIFO_H
#define DSY_SPECTRALANALYZERFIFO_H

#include <cstddef>

#include "spectral.h"
#include "shy_fft.h"


namespace daicsp
{
/** SpectralAnalyzerFifo
 * 
 *  Converts time-domain signals into frequency domain,
 *  and provides a FIFO interface for faster FFT 
 *  processing.
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-06-10
 *  Ported from Csound pvsanal.
 * 
 */
class SpectralAnalyzer
{
  public:
    SpectralAnalyzer() {}
    ~SpectralAnalyzer() {}

    /** Indicates the current status of the module. 
         *  Warnings are indicated by a leading W, and are silently corrected. 
         *  Errors are indicated by a leading E and cause an immediate exit.
         * 
         *  \param OK - No errors have been reported.
         *  \param E_FFT_NOT_POWER - Currently, the fftsize must be a power of 2.
         *  \param E_FFT_TOO_SMALL - The given fftsize is too small.
         *  \param W_OVERLAP_TOO_BIG - The frame overlap is too big (must be fftsize / 2 or less)
         *  \param E_OVERLAP_TOO_SMALL - The overlap must be greater than zero.
         *  \param E_BLOCK_TOO_BIG - The audio block in combination with the fft size is too large for the input buffer. 
         *  \param E_WINDOW_TOO_SMALL - The window must be equal to or greater than fftsize.
         *  \param W_INVALID_WINDOW - The window type is not valid (only Hamming and Hann are currently supported).
         *  \param W_BUFFER_UNDERFLOW - The input buffer was filled before the previous buffer was fully processed.
         *  \param W_INVALID_STATE - Against all odds, you've put the state_ property into an invalid state. 
         *  \param E_SLIDING_NOT_IMPLEMENTED - Sliding is currently not implemented, so the overlap size must be greater than the audio block size.
         */
    enum STATUS
    {
        OK = 0,
        E_FFT_NOT_POWER,
        E_FFT_TOO_SMALL,
        E_OVERLAP_TOO_SMALL,
        W_OVERLAP_TOO_BIG,
        E_BLOCK_TOO_BIG,
        E_WINDOW_TOO_SMALL,
        W_INVALID_WINDOW,
        W_BUFFER_UNDERFLOW,
        W_INVALID_STATE,
        E_SLIDING_NOT_IMPLEMENTED,
    };

    /** Initializes the SpectralAnalyzerFifo module.
         *  \param fft_size - This determines the size of the FFT. It must be a power of greater than 32.
         *  \param overlap_size - This determines the overlap between frequency frames. It should be at least fftsize / 4, but cannot be greater than fftsize / 2.
         *  \param window_size - This determines the size of the analysis window. It must be greater than or equal to fftsize.
         *  \param window_type - The windowing function. Currently, only Hamming and Hann are supported.
         *  \param sample_rate - The program sample rate.
         *  \param audio_block - The program audio block size.
         */
    void Init(uint32_t        fft_size,
              uint32_t        overlap_size,
              uint32_t        window_size,
              SPECTRAL_WINDOW window_type,
              size_t          sample_rate,
              size_t          audio_block); //pvsanalset

    /** Writes a single sample to the FIFO, and
     *  queues the bulk processing when appropriate.
     */
    void Process(float sample);

    /** Processes a block of incoming audio from the FIFO,
     *  blocks until the FIFO is filled.
     *  \returns - A reference to the internal `SpectralBuffer` containing the frequency-domain data.
     */
    SpectralBuffer& ParallelProcess(); // pvsanal

    SpectralBuffer& GetFsig() { return fsig_out_; }

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

    /** Returns the estimated latency in seconds
     * 
     */
    float GetEstimatedLatency() { return (float) (num_overlaps_ * overlap_) / sample_rate_; }

    /** Check if the analyzer is ready for processing
     * 
     */
    size_t GetQueuedFrameCount() { return input_count_; }

  private:
    /** Corresponds to pvsanal's pvssanalset -- Phase Vocoder Synthesis _sliding_ analysis set.
         *  This is not currently implemented, but can be useful for small overlap sizes.
         */
    void InitSliding(uint32_t        fft_size,
                     uint32_t        overlap_size,
                     uint32_t        window_size,
                     SPECTRAL_WINDOW window_type,
                     size_t          sample_rate,
                     size_t          block); // pvssanalset

    void ProcessSliding(const float* in, size_t size); // pvssanal

    void Tick(float sample); // anal_tick

    void GenerateFrame(); // generate_frame

    int    buflen_;
    float  RoverTwoPi_, TwoPioverR_, Fexact_;
    float* nextIn_;
    int    nI_, Ii_, IOi_;
    size_t    inptr_;

    float input_[kFFTMaxFrames];

    // This is how we manage the input FIFO, so it's larger than a single buffer.
    float overlapbuf_[kFFTMaxOverlapBuff];
    // This is the size of an individual overlap frame
    size_t overlap_;
    // This is the size of the rotating overlapbuf in terms of overlap frames
    size_t num_overlaps_;
    // This is the number of overlap frames that need to be processed
    size_t input_count_;
    float* input_segment_;
    float* process_segment_;

    float analbuf_[kFFTMaxFloats];
    float analbufOut_[kFFTMaxFloats];
    float analwinbuf_[kFFTMaxWindow];
    float oldInPhase_[kFFTMaxFloats];

    // If floats aren't enough quality, return to doubles
    // float trig_[FFT_SIZE];
    // float* cosine_;
    // float* sine_;

    // TODO -- return these to the above state for sliding
    float* trig_;
    float* cosine_;
    float* sine_;

    ShyFFT<float, kFFTMaxSize> fft_;

    SpectralBuffer fsig_out_;
    float          sample_rate_;
    STATUS         status_;
};

} // namespace daicsp

#endif // DSY_SPECTRALANALYZERFIFO_H
