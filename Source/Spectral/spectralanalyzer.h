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

        enum FFT {
            MAX_SIZE = 4096,
        };

        /** Initializes the SpectralAnalyzer module.
         *  \param p - description
         */
        void Init(float sampleRate);

        /** Processes a single sample and returns it.
         *  \param in - input sample
         */
        void Process(const float* in, size_t size, SpectralBuffer &fsig); // pvsanal (?)

        // TODO -- documentation and proper return types
        // NOTE -- these can be a private methods, right?
        void Analyze(size_t size); // pvssanal

        void Tick(float sample); // anal_tick

        void UpdateFrame(); // generate_frame

        void Interlace(float* fftSeparated, int length);


    private:
        // TODO -- convert these to proper types
        // OPDS    h; csound specific property that we don't need
        
        // This may change -- not necessary to be a member of this class
        SpectralBuffer  *fsig;          /* output signal is an analysis frame */
        float   *ain;                   /* input sig is audio */
        float   *fftsize;               /* params */
        float   *overlap;
        float   *winsize;
        float   *wintype;
        float   *format;                /* always PVS_AMP_FREQ at present */
        float   *init;                  /* not yet implemented */
        /* internal */
        int     buflen;
        float   fund,arate;
        float   RoverTwoPi,TwoPioverR,Fexact;
        float   *nextIn;
        int     nI,Ii,IOi;              /* need all these ?; double as N and NB */
        int     inptr;

        // TODO -- these can probably just be float arrays
        float*   input;
        float*   overlapbuf;
        float*   analbuf;
        float*   analbufOut;
        float*   analwinbuf;     /* prewin in SDFT case */
        float*   oldInPhase;
        float*   trig;
        float   *cosine, *sine; // If floats aren't enough quality, return to doubles
        // void    *setup; // this is a struct for Csound's FFT implementation
        ShyFFT<float, FFT::MAX_SIZE> fft_;

        float sr_;

        // TODO -- ideally, no swapping would be necessary
        float* swapBuffer_;

        // property to keep track of running input length
        int runningLength_;
        
};

} // namespace daicsp

#endif // DSY_SPECTRALANALYZER_H
