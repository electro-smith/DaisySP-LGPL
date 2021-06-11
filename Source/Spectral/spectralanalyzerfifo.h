#pragma once
#ifndef DSY_SPECTRALANALYZERFIFO_H
#define DSY_SPECTRALANALYZERFIFO_H

#include <cstddef>
#include <math.h>
#include <string.h>

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
 *  \param FFT_SIZE - This determines the size of the FFT. It must be a power of greater than 32.
 *  \param overlap - This determines the overlap between frequency frames. It should be at least fftsize / 4, but cannot be greater than fftsize / 2.
 *  \param windowSize - This determines the size of the analysis window. It must be greater than or equal to fftsize.
 */

template <size_t FFT_SIZE    = 2048,
          size_t OVERLAP     = 512,
          size_t WINDOW_SIZE = 2048>
class SpectralAnalyzerFifo
{
  public:
    SpectralAnalyzerFifo() {}
    ~SpectralAnalyzerFifo() {}

    /** Indicates the current status of the module. 
         *  Warnings are indicated by a leading W, and are silently corrected. 
         *  Errors are indicated by a leading E and cause an immediate exit.
         * 
         *  \param OK - No errors have been reported.
         *  \param E_FFT_NOT_POWER - Currently, the fftsize must be a power of 2.
         *  \param E_FFT_TOO_SMALL - The given fftsize is too small.
         *  \param W_OVERLAP_TOO_BIG - The frame overlap is too big (must be fftsize / 2 or less)
         *  \param E_OVERLAP_TOO_SMALL - The overlap must be greater than zero.
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
        E_WINDOW_TOO_SMALL,
        W_INVALID_WINDOW,
        W_BUFFER_UNDERFLOW,
        W_INVALID_STATE,
        E_SLIDING_NOT_IMPLEMENTED,
    };

    enum STATE
    {
        INIT = 0,
        IDLE,
        PROCESSING,
    };

    /** Initializes the SpectralAnalyzerFifo module.
         *  \param windowType - The windowing function. Currently, only Hamming and Hann are supported.
         *  \param sampleRate - The program sample rate.
         *  \param block - The program audio block size
         */
    void Init(SPECTRAL_WINDOW windowType,
              size_t          sampleRate,
              size_t          block); //pvsanalset

    /** Writes a single sample to the FIFO, and
     *  queues the bulk processing when appropriate.
     */
    void Sample(float sample);

    /** Processes a block of incoming audio from the FIFO,
     *  blocks until the FIFO is filled.
     *  \returns - A reference to the internal `SpectralBuffer` containing the frequency-domain data.
     */
    SpectralBuffer<FFT_SIZE + 2>& Process(); // pvsanal

    SpectralBuffer<FFT_SIZE + 2>& GetFsig() { return fsig_; }

    /** Retrieves the current status. Useful for error checking.
         */
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
    /** Corresponds to pvsanal's pvssanalset -- Phase Vocoder Synthesis _sliding_ analysis set.
         *  This is not currently implemented, but can be useful for small overlap sizes.
         */
    void InitSliding(SPECTRAL_WINDOW windowType,
                     size_t          sampleRate,
                     size_t          block); // pvssanalset

    void ProcessSliding(const float* in, size_t size); // pvssanal

    void Tick(float sample); // anal_tick

    void GenerateFrame(); // generate_frame

    /** Interlaces real and imaginary values from shy_fft.
         *  This is made necessary because Csound's fft function
         *  interleaves real and imaginary values by default.
         *  shy_fft, on the other hand, puts real in the first
         *  half and imaginary in the other.
         */
    void Interlace(float* fftSeparated, float* targetBuffer, const int length);

    int    buflen_;
    float  RoverTwoPi_, TwoPioverR_, Fexact_;
    float* nextIn_;
    int    nI_, Ii_, IOi_;
    int    inptr_;

    float input_[FFT_SIZE * 4];

    float overlapbuf_
        [OVERLAP
         * 2]; // This is how we manage the input FIFO, so it's twice the size
    size_t halfOverlap_;
    float* inputSegment_;
    float* processSegment_;

    float analbuf_[FFT_SIZE + 2];
    float analbufOut_[FFT_SIZE + 2];
    float analwinbuf_[FFT_SIZE + 1];
    float oldInPhase_[FFT_SIZE / 2 + 1];

    // If floats aren't enough quality, return to doubles
    // float trig_[FFT_SIZE];
    // float* cosine_;
    // float* sine_;

    // TODO -- return these to the above state for sliding
    float* trig_;
    float* cosine_;
    float* sine_;

    ShyFFT<float, FFT_SIZE> fft_;

    SpectralBuffer<FFT_SIZE + 2> fsig_;
    float                        sr_;
    STATUS                       status_;
    STATE                        state_;
};

#include "spectralanalyzerfifoimpl.h"

} // namespace daicsp

#endif // DSY_SPECTRALANALYZERFIFO_H
