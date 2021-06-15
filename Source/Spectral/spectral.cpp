#include <cmath>
#include "dsp.h"
#include "spectral.h"

using namespace daicsp;
using namespace std;

int daicsp::GetPasses(int fft_size)
{
    int fft_num_passes_ = 0;
    for(size_t t = fft_size; t > 1; t >>= 1)
    {
        ++fft_num_passes_;
    }
    return fft_num_passes_;
}

float daicsp::mod2Pi(float value)
{
    value = fmod(value, TWOPI_F);
    if(value <= -PI_F)
    {
        return value + TWOPI_F;
    }
    else if(value > PI_F)
    {
        return value - TWOPI_F;
    }
    else
        return value;
}

void daicsp::SpectralWindow(float* windowBuffer, int type, int windowLength)
{
    // float  fpos, inc;
    // float   *ftable;
    // int     i, n, flen, even;
    int i, n, even;

    even = (windowLength + 1) & 1;
    switch(type)
    {
        case SPECTRAL_WINDOW::HAMMING:
            Hamming(windowBuffer, (windowLength >> 1), even);
            // return OK;
            return;
        case SPECTRAL_WINDOW::HANN: /* Hanning */
            Vonhann(windowBuffer, (windowLength >> 1), even);
            // return OK;
            return;
        case SPECTRAL_WINDOW::KAISER: /* Kaiser */
        {
            float beta = 6.8f;
            float x, flen2, besbeta;
            flen2 = 1.0f
                    / ((float)(windowLength >> 1) * (float)(windowLength >> 1));
            besbeta = 1.0f / Besseli(beta);
            n       = windowLength >> 1;
            x       = (even ? 0.5f : 0.05f);
            for(i = 0; i < n; i++, x += 1.0f)
                windowBuffer[i]
                    = (float)(Besseli(beta * sqrt(1.0f - x * x * flen2))
                              * besbeta);
            windowBuffer[i] = 0.0f;
        }
            // return OK;
            return;
        default:
            // TODO -- error handling
            if(type >= 0)
            {
                // return csound->InitError(csound, Str("invalid window type"));
                return;
            }
    }

    // TODO -- permit user-generated windows?

    // /* use table created with GEN20 */
    // flen = csoundGetTable(csound, &ftable, -(type));
    // if (UNLIKELY(flen < 0))
    //   return csound->InitError(csound, Str("ftable for window not found"));
    // inc = (float)flen / (float)(windowLength & (~1));
    // fpos = ((float)flen + (float)even * inc) * 0.5;
    // n = windowLength >> 1;
    // /* this assumes that for a window with even size, space for an extra */
    // /* sample is allocated */
    // for (i = 0; i < n; i++) {
    //   float  frac, tmp;
    //   int     pos;
    //   frac = modf(fpos, &tmp);
    //   pos = (int) tmp;
    //   windowBuffer[i] = ftable[pos] + ((ftable[pos + 1] - ftable[pos]) * (float) frac);
    //   fpos += inc;
    // }
    // windowBuffer[n] = (even ? FL(0.0) : ftable[flen]);
    // // return OK;
}

void daicsp::Hamming(float* windowBuffer, int windowLength, int even)
{
    float ftmp;
    int   i;

    ftmp = PI_F / windowLength;

    if(even)
    {
        for(i = 0; i < windowLength; i++)
            windowBuffer[i] = (0.54f + 0.46f * cos(ftmp * ((float)i + 0.5f)));
        windowBuffer[windowLength] = 0.0f;
    }
    else
    {
        windowBuffer[0] = 1.0f;
        for(i = 1; i <= windowLength; i++)
            windowBuffer[i] = (0.54f + 0.46f * cos(ftmp * (float)i));
    }
}

void daicsp::Vonhann(float* windowBuffer, int windowLength, int even)
{
    float ftmp;
    int   i;

    ftmp = PI_F / windowLength;

    if(even)
    {
        for(i = 0; i < windowLength; i++)
            windowBuffer[i]
                = (float)(0.5 + 0.5 * cos(ftmp * ((float)i + 0.5f)));
        windowBuffer[windowLength] = 0.0f;
    }
    else
    {
        windowBuffer[0] = 1.0f;
        for(i = 1; i <= windowLength; i++)
            windowBuffer[i] = (float)(0.5 + 0.5 * cos(ftmp * (float)i));
    }
}

float daicsp::Besseli(float x)
{
    float ax, ans;
    float y;

    if((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        ans = (1.0
               + y
                     * (3.5156229f
                        + y
                              * (3.0899424f
                                 + y
                                       * (1.2067492f
                                          + y
                                                * (0.2659732f
                                                   + y
                                                         * (0.360768e-1f
                                                            + y * 0.45813e-2f))))));
    }
    else
    {
        y   = 3.75 / ax;
        ans = ((exp(ax) / sqrt(ax))
               * (0.39894228f
                  + y
                        * (0.1328592e-1f
                           + y
                                 * (0.225319e-2f
                                    + y
                                          * (-0.157565e-2f
                                             + y
                                                   * (0.916281e-2f
                                                      + y
                                                            * (-0.2057706e-1f
                                                               + y
                                                                     * (0.2635537e-1f
                                                                        + y
                                                                              * (-0.1647633e-1f
                                                                                 + y * 0.392377e-2f)))))))));
    }
    return ans;
}

void daicsp::Interleave(float*    fft_separated,
                        float*    target_buffer,
                        const int length)
{
    // unfortunately, interleaving in place is not trivial, so another buffer will have to do
    int halflen = length / 2;
    for(int i = 0; i < halflen; i++)
    {
        target_buffer[i * 2]     = fft_separated[i];
        target_buffer[i * 2 + 1] = fft_separated[i + halflen];
    }
}

void daicsp::Deinterleave(float*    fft_interlaced,
                          float*    working_buffer,
                          const int length)
{
    int halflen = length / 2;

    for(int i = 0; i < halflen; i++)
    {
        working_buffer[i]           = fft_interlaced[i * 2];
        working_buffer[i + halflen] = fft_interlaced[i * 2 + 1];
    }

    for(int i = 0; i < length; i++)
    {
        fft_interlaced[i] = working_buffer[i];
    }
}