#pragma once
#ifndef DSY_SPECTRAL_H
#define DSY_SPECTRAL_H

namespace daicsp
{

enum WINDOW_SIZE {
    MAX = 4096,
};

enum SPECTRAL_WINDOW {
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

enum SPECTRAL_FORMAT {
    AMP_FREQ = 0,
    AMP_PHASE,
    COMPLEX,
    TRACKS,
};

typedef struct {
    int N;
    bool sliding;
    int NB;
    int overlap;
    int winsize;
    int wintype;
    int format;
    unsigned int framecount;
    float frame[WINDOW_SIZE::MAX];

    void Init(SPECTRAL_WINDOW window, int windowSize);

} SpectralBuffer;


typedef struct {
    float a;
    float b;
} Complex;


float mod2Pi(float value);


} // namespace daicsp

#endif