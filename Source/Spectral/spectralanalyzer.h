#pragma once
#ifndef DSY_SPECTRALANALYZER_H
#define DSY_SPECTRALANALYZER_H

#include "spectral.h"
#include "shy_fft.h"

namespace daicsp
{

/** SpectralAnalyzer
 *  Author: Gabriel Ball
 *  Date: 2021-06-03
 *  
 *  Ported from Csound pvs opcodes.
 */


class SpectralAnalyzer
{
    public:
        SpectralAnalyzer () {}
        ~SpectralAnalyzer () {}

        /** Initializes the SpectralAnalyzer module.
         *  \param p - description
         */
        void Init(int fftsize, int overlap, int windowSize, SPECTRAL_WINDOW windowType, size_t sampleRate, size_t block); //pvsanalset

        /** Processes a single sample and returns it.
         *  \param in - input sample
         */
        SpectralBuffer& Process(const float* in, size_t size); // pvsanal

        SpectralBuffer& GetFsig() { return fsig_; }
        
    private:

        void InitSmall(); // pvssanalset

        // TODO -- documentation and proper return types
        // NOTE -- these can be a private methods, right?
        void Analyze(const float* in, size_t size); // pvssanal

        void Tick(float sample); // anal_tick

        void UpdateFrame(); // generate_frame

        void Interlace(float* fftSeparated, float* targetBuffer, const int length);


    private:
        // TODO -- convert these to proper types
        // OPDS    h; csound specific property that we don't need
        
        // This may change -- not necessary to be a member of this class
        // SpectralBuffer  *fsig;          /* output signal is an analysis frame */
        // float   *ain;                   /* input sig is audio */
        // float   *fftsize;               /* params */
        // float   *overlap;
        // float   *winsize;
        // float   *wintype;
        // float   *format;                /* always PVS_AMP_FREQ at present */
        // float   *init;                  /* not yet implemented */
        /* internal */
        int     buflen;
        // float   fund,arate;
        float   RoverTwoPi,TwoPioverR,Fexact;
        float   *nextIn;
        int     nI,Ii,IOi;              /* need all these ?; double as N and NB */
        int     inptr;

        float input[WINDOW_SIZE::MAX];
        float overlapbuf[WINDOW_SIZE::MAX];
        float analbuf[WINDOW_SIZE::MAX];
        float analbufOut[WINDOW_SIZE::MAX];
        float analwinbuf[WINDOW_SIZE::MAX];     /* prewin in SDFT case */
        float oldInPhase[WINDOW_SIZE::MAX];

        // TODO -- adjust these according to their appropriate sizes
        float trig[WINDOW_SIZE::MAX];
        float cosine[WINDOW_SIZE::MAX];
        float sine[WINDOW_SIZE::MAX]; // If floats aren't enough quality, return to doubles

        ShyFFT<float, WINDOW_SIZE::MAX> fft_;
        SpectralBuffer fsig_;
        float sr_;

        // // TODO -- ideally, no swapping would be necessary
        // float* swapBuffer_;
        
};

} // namespace daicsp

#endif // DSY_SPECTRALANALYZER_H
