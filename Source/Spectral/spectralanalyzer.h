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
        void Init(float sampleRate);

        /** Processes a single sample and returns it.
         *  \param in - input sample
         */
        void Process(const float* in, size_t size, SpectralBuffer &fsig); // pvsanal (?)

        // TODO -- documentation and proper return types
        // NOTE -- these can be a private methods, right?
        void Analyze(const float* in, size_t size, SpectralBuffer &fsig); // pvssanal

        void Tick(float sample, SpectralBuffer &fsig); // anal_tick

        void UpdateFrame(SpectralBuffer &fsig); // generate_frame

        void Interlace(float* fftSeparated, int length);


    private:
        // TODO -- convert these to proper types
        // OPDS    h; csound specific property that we don't need
        
        // This may change -- not necessary to be a member of this class
        // SpectralBuffer  *fsig;          /* output signal is an analysis frame */
        // float   *ain;                   /* input sig is audio */
        // float   *fftsize;               /* params */
        float   *overlap;
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

        // TODO -- these can probably just be float arrays
        // TODO -- how big should they be, though??
        float input[WINDOW_SIZE::MAX];
        float overlapbuf[WINDOW_SIZE::MAX];
        float analbuf[WINDOW_SIZE::MAX];
        float analbufOut[WINDOW_SIZE::MAX];
        float analwinbuf[WINDOW_SIZE::MAX];     /* prewin in SDFT case */
        float oldInPhase[WINDOW_SIZE::MAX];
        float trig[WINDOW_SIZE::MAX];
        float cosine[WINDOW_SIZE::MAX];
        float sine[WINDOW_SIZE::MAX]; // If floats aren't enough quality, return to doubles
        // void    *setup; // this is a struct for Csound's FFT implementation
        ShyFFT<float, WINDOW_SIZE::MAX> fft_;

        float sr_;

        // TODO -- ideally, no swapping would be necessary
        float* swapBuffer_;

        // property to keep track of running input length
        int runningLength_;
        
};

} // namespace daicsp

#endif // DSY_SPECTRALANALYZER_H
