#pragma once
#ifndef DSY_SPECTRALANALYZER_H
#define DSY_SPECTRALANALYZER_H

#include "spectral.h"

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
        void Init();

        /** Processes a single sample and returns it.
         *  \param in - input sample
         */
        float Process(const float &in); // pvsanal (?)

        // TODO -- documentation and proper return types
        // NOTE -- these can be a private methods, right?
        void Analyze(float sample); // pvssanal

        void Tick(float sample); // anal_tick

        void UpdateBuffer(); // generate_frame

        float mod2Pi(float value); // mod2Pi


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
        AUXCH   input;
        AUXCH   overlapbuf;
        AUXCH   analbuf;
        AUXCH   analwinbuf;     /* prewin in SDFT case */
        AUXCH   oldInPhase;
        AUXCH   trig;
        float   *cosine, *sine; // If floats aren't enough quality, return to doubles
        void    *setup; // this is a struct for Csound's FFT implementation
        
};

} // namespace daicsp

#endif // DSY_SPECTRALANALYZER_H
