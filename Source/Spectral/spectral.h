#pragma once
#ifndef DSY_SPECTRAL_H
#define DSY_SPECTRAL_H

namespace daicsp
{
#define __MAX_FFT__ 4096

enum FFT
{
    MAX_SIZE    = __MAX_FFT__,
    MAX_WINDOW  = __MAX_FFT__ + 1,
    MAX_FLOATS  = __MAX_FFT__ + 2,
    MAX_BINS    = __MAX_FFT__ / 2 + 1,
    MAX_OVERLAP = __MAX_FFT__ / 2,
    MAX_FRAMES  = __MAX_FFT__ * 4,

    // NOTE -- This is a temporary measure to ensure
    // the appropriate buffers have enough space for
    // sliding and non-sliding applications. It is
    // derived from: numBins * audioBlock
    // Needless to say, this is a temporary hack and
    // needs a better solution since the audio block
    // size may not be known at compile time.
    // (this assumes 48)
    MAX_SLIDE_SIZE = 4704,
};

enum SPECTRAL_WINDOW
{
    HAMMING = 0,
    HANN,
    KAISER,
    CUSTOM,
    BLACKMAN,
    BLACKMAN_EXACT,
    NUTTALLC3,
    BHARRIS_3,
    BHARRIS_MIN,
    RECT,
};

enum SPECTRAL_FORMAT
{
    AMP_FREQ = 0,
    AMP_PHASE,
    COMPLEX,
    TRACKS,
};

typedef struct
{
    int          N;
    bool         sliding;
    int          NB;
    int          overlap;
    int          winsize;
    int          wintype;
    int          format;
    unsigned int framecount;
    float        frame[FFT::MAX_FLOATS];
} SpectralBuffer;


typedef struct
{
    float a;
    float b;
} Complex;

int GetPasses(int fft_size);

float mod2Pi(float value);

void SpectralWindow(float* windowBuffer, int windowType, int windowLength);

void Hamming(float* windowBuffer, int windowLength, int even);

void Vonhann(float* windowBuffer, int windowLength, int even);

double Besseli(double x);


} // namespace daicsp

#endif