#pragma once
#ifndef DSY_SPECTRAL_H
#define DSY_SPECTRAL_H

#include <cstddef>

namespace daicsp
{
// NOTE -- this macro determines the static memory allocated
// for all spectral arrays. If you'd like to minimize the
// memory footprint, simply define `__FFT_SIZE__` before including
// daicsp
#ifndef __FFT_SIZE__
#define __FFT_SIZE__ 4096
#endif

static constexpr unsigned int kFFTMaxSize    = __FFT_SIZE__;
static constexpr unsigned int kFFTMaxWindow  = __FFT_SIZE__ + 1;
static constexpr unsigned int kFFTMaxFloats  = __FFT_SIZE__ + 2;
static constexpr unsigned int kFFTMaxBins    = __FFT_SIZE__ / 2 + 1;
static constexpr unsigned int kFFTMaxOverlap = __FFT_SIZE__ / 2;
static constexpr unsigned int kFFTMaxFrames  = __FFT_SIZE__ * 4;
// NOTE -- This is a temporary measure to ensure
// the appropriate buffers have enough space for
// sliding and non-sliding applications. It is
// derived from: numBins * audioBlock
// Needless to say, this is a temporary hack and
// needs a better solution since the audio block
// size may not be known at compile time.
// (this assumes 48)
static constexpr unsigned int kFFTMaxSlide = 4704;

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

struct SpectralBuffer
{
    int          N;
    bool         sliding;
    int          NB;
    int          overlap;
    int          winsize;
    int          wintype;
    int          format;
    unsigned int framecount;
    float        frame[kFFTMaxFloats];
    bool         ready;
};

// This causes errors somehow
// Allowing for easy test to ensure fsigIn != fsigOut
// bool operator==(SpectralBuffer& fsig1, SpectralBuffer& fsig2)
// {
//     if(&fsig1 == &fsig2) return true;
//     else return false;
// }
// bool operator!=(SpectralBuffer& fsig1, SpectralBuffer& fsig2)
// {
//     if(&fsig1 != &fsig2) return true;
//     else return false;
// }

typedef struct
{
    float real;
    float imaginary;
} Complex;

int GetPasses(int fft_size);

float mod2Pi(float value);

void SpectralWindow(float* windowBuffer, int windowType, int windowLength);

void Hamming(float* windowBuffer, int windowLength, int even);

void Vonhann(float* windowBuffer, int windowLength, int even);

float Besseli(float x);

/** Interlaces real and imaginary values from shy_fft.
         *  This is made necessary because Csound's fft function
         *  interleaves real and imaginary values by default.
         *  shy_fft, on the other hand, puts real in the first
         *  half and imaginary in the other.
         */
void Interleave(float* fft_separated, float* target_buffer, const int length);

/** Deinterlaces real and imaginary values into blocks for processing by shy_fft.
         */
void Deinterleave(float* interlaced, float* target_buffer, const int length);


} // namespace daicsp

#endif