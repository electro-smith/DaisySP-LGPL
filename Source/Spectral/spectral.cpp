#include <math.h>
#include "dsp.h"
#include "spectral.h"

using namespace daicsp;

float daicsp::mod2Pi(float value) {
    value = fmod(value,TWOPI_F);
    if (value <= -PI_F) {
        return value + TWOPI_F;
    }
    else if (value > PI_F) {
        return value - TWOPI_F;
    }
    else
      return value;
}

void SpectralBuffer::Init(SPECTRAL_WINDOW window, int windowSize) {
    // TODO -- not sure how these should be set
    sliding = false;
    overlap = 0;
    framecount = 0;
    N = windowSize; // almost certainly wrong
    NB = N;

    // These are probably right
    wintype = (int) window;
    winsize = windowSize;
    // NOTE -- the PVANAL struct mentions that the format is always AMP_FREQ
    // at present, so we'll set that here by default;
    format = (int) SPECTRAL_FORMAT::AMP_FREQ;
}