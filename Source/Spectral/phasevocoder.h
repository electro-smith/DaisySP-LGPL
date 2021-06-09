#pragma once
#ifndef DSY_PHASEVOCODER_H
#define DSY_PHASEVOCODER_H

#include "spectral.h"
#include "shy_fft.h"

namespace daicsp
{

/** PhaseVocoder
 * 
 *  Converts frequency-domain signals back to time-domain.
 * 
 *  Author: Gabriel Ball
 *  Date: 2021-06-03
 */
class PhaseVocoder
{
    public:
        PhaseVocoder () {}
        ~PhaseVocoder () {}

        enum STATUS {
            OK = 0,
            SLIDING_NOT_IMPLEMENTED,

        };
       
        // NOTE -- this might be useful as an overload?
        // void Init(int fftsize, int overlap, int windowSize, SPECTRAL_WINDOW windowType, size_t sampleRate, size_t block);

        /** Initializes the PhaseVocoder module.
         *  \param fsig - Initialized frequency-domain signal from the intended source.
         *  \param sampleRate - The program sample rate.
         *  \param block - The program audio block size.
         */
        void Init(SpectralBuffer fsig, size_t sampleRate, size_t block);

        // NOTE -- probably best to not return a raw pointer
        /** Processes an incoming fsig and returns a pointer to the resynthesized audio of the given size.
         *  \param fsig - A `SpectralBuffer` from the source used to initialize this module.
         *  \param size - The size of the audio block.
         *  \returns - Pointer to an internal buffer of resynthesizes audio
         */
        float* Process(SpectralBuffer &fsig, size_t size); // pvsynth

        /** Retrieves the current status. Useful for error checking.
         */
        STATUS GetStatus() { return status_; } 

    private:

        /** Corresponds to pvsynth's pvssynthset -- Phase Vocoder Synthesis _sliding_ set.
         *  This is not currently implemented, but can be useful for small overlap sizes.
         */
        void ProcessSliding(SpectralBuffer& fsig, size_t size); // pvssynth

        float Tick(SpectralBuffer &fsig); // analyze_tick

        void GenerateFrame(SpectralBuffer &fsig); // process_frame

        /** Interlaces real and imaginary values into real and imaginary blocks for use with shy_fft.
         *  This is made necessary because Csound's fft function
         *  Interleaves real and imaginary values by default.
         *  shy_fft, on the other hand, puts real in the first
         *  half and imaginary in the other.
         */
        void Deinterlace(float* interlaced, float* targetBuffer, const int length);

        int     buflen_;
        float   RoverTwoPi_,TwoPioverR_,Fexact_;
        float   *nextOut_;
        int     nO_,Ii_,IOi_;
        int     outptr_;

        float output_[FFT::MAX_FRAMES];
        float overlapbuf_[FFT::MAX_OVERLAP];
        float synbuf_[FFT::MAX_FLOATS];
        float synbufOut_[FFT::MAX_FLOATS];
        float analwinbuf_[FFT::MAX_WINDOW]; 
        float synwinbuf_[FFT::MAX_WINDOW];
        float oldOutPhase_[FFT::MAX_BINS];

        ShyFFT<float, FFT::MAX_SIZE> fft_;
        float sr_;
        int blockSize_;

        // NOTE -- ensure this is always greater than the block size!
        float outputBuffer_[64];
        /* check these against fsig vals */
        // int overlap,winsize,fftsize,wintype,format;
        int bin_index_;      /* for phase normalization across frames */
        STATUS status_;
};

} // namespace daicsp

#endif // DSY_PHASEVOCODER_H
