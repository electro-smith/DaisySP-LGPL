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
          size_t OVERLAP     = 512,
          size_t WINDOW_SIZE = 2048>
class PhaseVocoderFifo
{
  public:
    PhaseVocoderFifo() {}
    ~PhaseVocoderFifo() {}

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
         *  \param fsig - Initialized frequency-domain signal from the intended source.
         *  \param sampleRate - The program sample rate.
         *  \param block - The program audio block size.
         */
    void
    Init(SpectralBuffer<FFT_SIZE + 2>& fsig, size_t sampleRate, size_t block);

    // NOTE -- probably best to not return a raw pointer
    /** Processes an incoming fsig and returns a pointer to the resynthesized audio of the given size.
         *  \param fsig - A `SpectralBuffer` from the source used to initialize this module.
         *  \param size - The size of the audio block.
         *  \returns - Pointer to an internal buffer of resynthesizes audio
         */
    void Process(SpectralBuffer<FFT_SIZE + 2>& fsig); // pvsynth

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
            while(1)
                ;
    }

  private:
    /** Corresponds to pvsynth's pvssynthset -- Phase Vocoder Synthesis _sliding_ set.
         *  This is not currently implemented, but can be useful for small overlap sizes.
         */
    void ProcessSliding(SpectralBuffer<FFT_SIZE + 2>& fsig,
                        size_t                        size); // pvssynth

    float Tick(SpectralBuffer<FFT_SIZE + 2>& fsig); // analyze_tick

    void GenerateFrame(SpectralBuffer<FFT_SIZE + 2>& fsig); // process_frame

    /** Interlaces real and imaginary values into real and imaginary blocks for use with shy_fft.
         *  This is made necessary because Csound's fft function
         *  Interleaves real and imaginary values by default.
         *  shy_fft, on the other hand, puts real in the first
         *  half and imaginary in the other.
         */
    void Deinterlace(float* interlaced, float* targetBuffer, const int length);

    int    buflen_;
    float  RoverTwoPi_, TwoPioverR_, Fexact_;
    float* nextOut_;
    int    nO_, Ii_, IOi_;
    size_t outptr_;

    float  output_[FFT_SIZE * 4];
    float  overlapbuf_[OVERLAP * 2];
    size_t halfOverlap_;
    float* outputSegment_;
    float* processSegment_;

    float synbuf_[FFT_SIZE + 2];
    float synbufOut_[FFT_SIZE + 2];
    float analwinbuf_[FFT_SIZE + 1];
    float synwinbuf_[FFT_SIZE + 1];
    float oldOutPhase_[FFT_SIZE / 2 + 1];

    ShyFFT<float, FFT_SIZE> fft_;

    float sr_;
    int   blockSize_;

    /* check these against fsig vals */
    // int overlap,winsize,fftsize,wintype,format;
    int    bin_index_; /* for phase normalization across frames */
    STATUS status_;
    STATE  state_;
};

#include "phasevocoderfifoimpl.h"

} // namespace daicsp

#endif // DSY_PHASEVOCODERFIFO_H
