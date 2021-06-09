#pragma once
#ifndef DSY_SPECTRALANALYZER_H
#define DSY_SPECTRALANALYZER_H

#include "spectral.h"
#include "shy_fft.h"

namespace daicsp
{

/** SpectralAnalyzer
 * 
 *  Converts time-domain signals into frequency domain.
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-06-03
 *  Ported from Csound pvsanal.
 */


class SpectralAnalyzer
{
    public:
        SpectralAnalyzer () {}
        ~SpectralAnalyzer () {}

        /** Indicates the current status of the module. 
         *  Warnings are indicated by a leading W, and are silently corrected. 
         *  Errors are indicated by a leading E and cause an immediate exit.
         * 
         *  \param OK - No errors have been reported.
         *  \param E_FFT_NOT_POWER - Currently, the fftsize must be a power of 2.
         *  \param W_FFT_TOO_SMALL - The given fftsize is too small.
         *  \param W_FFT_TOO_BIG - The given fftsize is too big, limited by memory usage.
         *  \param W_OVERLAP_TOO_BIG - The frame overlap is too big (must be fftsize / 2 or less)
         *  \param E_OVERLAP_TOO_SMALL - The overlap must be greater than zero.
         *  \param W_WINDOW_TOO_SMALL - The window must be equal to or greater than fftsize.
         *  \param W_WINDOW_TOO_BIG - The window is too big, limited by memory usage.
         *  \param W_INVALID_WINDOW - The window type is not valid (only Hamming and Hann are currently supported).
         *  \param E_SLIDING_NOT_IMPLEMENTED - Sliding is currently not implemented, so the overlap size must be greater than the audio block size.
         */
        enum STATUS {
            OK = 0,
            E_FFT_NOT_POWER,
            W_FFT_TOO_SMALL,
            W_FFT_TOO_BIG,
            E_OVERLAP_TOO_SMALL,
            W_OVERLAP_TOO_BIG,
            W_WINDOW_TOO_SMALL,
            W_WINDOW_TOO_BIG,
            W_INVALID_WINDOW,
            E_SLIDING_NOT_IMPLEMENTED,
        };

        /** Initializes the SpectralAnalyzer module.
         *  \param fftsize - This determines the size of the FFT. It must be a power of two between 64 and 2048 inclusive.
         *  \param overlap - This determines the overlap between frequency frames. It should be at least fftsize / 4, but cannot be greater than fftsize / 2.
         *  \param windowSize - This determines the size of the analysis window. It must be greater than or equal to fftsize.
         *  \param windowType - The windowing function. Currently, only Hamming and Hann are supported.
         *  \param sampleRate - The program sample rate.
         *  \param block - The program audio block size
         */
        void Init(int fftsize, int overlap, int windowSize, SPECTRAL_WINDOW windowType, size_t sampleRate, size_t block); //pvsanalset

        /** Processes a block of incoming audio.
         *  \param in - Pointer to a buffer of audio data.
         *  \param size - The size of the given buffer.
         *  \returns - A reference to the internal `SpectralBuffer` containing the frequency-domain data.
         */
        SpectralBuffer& Process(const float* in, size_t size); // pvsanal

        SpectralBuffer& GetFsig() { return fsig_; }

        /** Retrieves the current status. Useful for error checking.
         */
        STATUS GetStatus() { return status_; } 
        
    private:

        /** Corresponds to pvsanal's pvssanalset -- Phase Vocoder Synthesis _sliding_ analysis set.
         *  This is not currently implemented, but can be useful for small overlap sizes.
         */
        void InitSliding(int fftsize, int overlap, int windowSize, SPECTRAL_WINDOW windowType, size_t sampleRate, size_t block); // pvssanalset

        void ProcessSliding(const float* in, size_t size); // pvssanal

        void Tick(float sample); // anal_tick

        void GenerateFrame(); // generate_frame

        /** Interlaces real and imaginary values from shy_fft.
         *  This is made necessary because Csound's fft function
         *  Interleaves real and imaginary values by default.
         *  shy_fft, on the other hand, puts real in the first
         *  half and imaginary in the other.
         */
        void Interlace(float* fftSeparated, float* targetBuffer, const int length);

        int     buflen_;
        float   RoverTwoPi_,TwoPioverR_,Fexact_;
        float   *nextIn_;
        int     nI_,Ii_,IOi_;  
        int     inptr_;

        float input_[FFT::MAX_FRAMES];
        float overlapbuf_[FFT::MAX_OVERLAP];
        float analbuf_[FFT::MAX_FLOATS];
        float analbufOut_[FFT::MAX_FLOATS];
        float analwinbuf_[FFT::MAX_WINDOW];
        float oldInPhase_[FFT::MAX_BINS];

        // If floats aren't enough quality, return to doubles
        // float trig_[FFT::MAX_SIZE];
        // float* cosine_;
        // float* sine_; 

        // TODO -- return these to the above state for sliding
        float* trig_;
        float* cosine_;
        float* sine_;

        ShyFFT<float, FFT::MAX_SIZE> fft_;
        SpectralBuffer fsig_;
        float sr_;
        STATUS status_;
};

} // namespace daicsp

#endif // DSY_SPECTRALANALYZER_H
