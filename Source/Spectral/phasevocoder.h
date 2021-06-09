#pragma once
#ifndef DSY_PHASEVOCODER_H
#define DSY_PHASEVOCODER_H

#include "spectral.h"
#include "shy_fft.h"

namespace daicsp
{

/** PhaseVocoder
 *  Author: Gabriel Ball
 *  Date: 2021-06-03
 */
class PhaseVocoder
{
    public:
        PhaseVocoder () {}
        ~PhaseVocoder () {}

        /** Initializes the PhaseVocoder module.
         *  \param p - description
         */
        // TODO -- this might be useful as an overload?
        // void Init(int fftsize, int overlap, int windowSize, SPECTRAL_WINDOW windowType, size_t sampleRate, size_t block);

        void Init(SpectralBuffer fsig, size_t sampleRate, size_t block);

        /** Processes a single sample and returns it.
         *  \param in - input sample
         */
        // TODO -- probably best to not pass a raw pointer
        float* Process(SpectralBuffer &fsig, size_t size); // pvsynth

        void Analyze(SpectralBuffer& fsig, size_t size); // pvssynth

        float Tick(SpectralBuffer &fsig); // analyze_tick

        void UpdateFrame(SpectralBuffer &fsig); // process_frame

        void Deinterlace(float* interlaced, float* targetBuffer, const int length);

    private:

        int     buflen;
        // float   fund,arate;
        float   RoverTwoPi,TwoPioverR,Fexact;
        float   *nextOut;
        int     nO,Ii,IOi;              /* need all these ?; double as N and NB */
        int     outptr;

        // TODO -- adjust these according to the MAX window size
        float output[WINDOW_SIZE::MAX];
        float overlapbuf[WINDOW_SIZE::MAX];
        float synbuf[WINDOW_SIZE::MAX];
        float synbufOut[WINDOW_SIZE::MAX];
        float analwinbuf[WINDOW_SIZE::MAX];     /* prewin in SDFT case */
        float synwinbuf[WINDOW_SIZE::MAX];
        float oldOutPhase[WINDOW_SIZE::MAX];

        ShyFFT<float, WINDOW_SIZE::MAX> fft_;
        float sr_;
        int blockSize_;

        // TODO -- ensure this is always greater than the block size!
        float outputBuffer[64];
        /* check these against fsig vals */
        // int overlap,winsize,fftsize,wintype,format;
        int bin_index;      /* for phase normalization across frames */
};

} // namespace daicsp

#endif // DSY_PHASEVOCODER_H
