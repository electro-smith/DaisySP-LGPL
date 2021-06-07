#pragma once
#ifndef DSY_PHASEVOCODER_H
#define DSY_PHASEVOCODER_H

#include "spectral.h"

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
        void Init(float sampleRate);

        /** Processes a single sample and returns it.
         *  \param in - input sample
         */
        void Process(const float* in, size_t size, SpectralBuffer &fsig); // pvsynth

        void Analyze(float sample); // pvssynth

        void Tick(float sample, SpectralBuffer &fsig); // analyze_tick

        void UpdateFrame(SpectralBuffer &fsig); // process_frame


    private:
        // OPDS    h;
        float   *aout;                  /* audio output signal */
        SpectralBuffer  *fsig;                  /* input signal is an analysis frame */
        float   *init;                  /* not yet implemented */
        /* internal */
        /* check these against fsig vals */
        int    overlap,winsize,fftsize,wintype,format;
        /* can we allow variant window tpes?  */
        int    buflen;
        float   fund,arate;
        float   RoverTwoPi,TwoPioverR,Fexact;
        float   *nextOut;
        int    nO,Ii,IOi;      /* need all these ?*/
        int    outptr;
        int    bin_index;      /* for phase normalization across frames */
        /* renderer gets all format info from fsig */

        float*   output;
        float*   overlapbuf;
        float*   synbuf;
        float*   analwinbuf;     /* may get away with a local alloc and free */
        float*   synwinbuf;
        float*   oldOutPhase;

        void    *setup; // FFT setup pointer

        // TODO -- private properties need trailing underscores
        float sr_;
};

} // namespace daicsp

#endif // DSY_PHASEVOCODER_H
