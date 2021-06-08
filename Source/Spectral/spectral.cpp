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

void daicsp::SpectralWindow(float* windowBuffer, int type, int windowLength) {
    // double  fpos, inc;
    // float   *ftable;
    // int     i, n, flen, even;
    int i, n, even;

    even = (windowLength + 1) & 1;
    switch (type) {
      case SPECTRAL_WINDOW::HAMMING:
        Hamming(windowBuffer, (windowLength >> 1), even);
        // return OK;
        return;
      case SPECTRAL_WINDOW::HANN:   /* Hanning */
        Vonhann(windowBuffer, (windowLength >> 1), even);
        // return OK;
        return;
      case SPECTRAL_WINDOW::KAISER:   /* Kaiser */
        {
          double  beta = 6.8;
          double  x, flen2, besbeta;
          flen2 = 1.0 / ((double)(windowLength >> 1) * (double)(windowLength >> 1));
          besbeta = 1.0 / Besseli(beta);
          n = windowLength >> 1;
          x = (even ? 0.5 : 0.05);
          for (i = 0; i < n; i++, x += 1.0)
            windowBuffer[i] = (float)(Besseli(beta * sqrt(1.0 - x * x * flen2))
                              * besbeta);
          windowBuffer[i] = 0.0f;
        }
        // return OK;
        return;
      default:
        // TODO -- error handling
        if (type >= 0) 
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
    // inc = (double)flen / (double)(windowLength & (~1));
    // fpos = ((double)flen + (double)even * inc) * 0.5;
    // n = windowLength >> 1;
    // /* this assumes that for a window with even size, space for an extra */
    // /* sample is allocated */
    // for (i = 0; i < n; i++) {
    //   double  frac, tmp;
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
    double ftmp;
    int i;

    ftmp = PI_F/windowLength;

    if (even) {
      for (i=0; i<windowLength; i++)
        windowBuffer[i] = (float)(0.54 + 0.46*cos(ftmp*((double)i+0.5)));
      windowBuffer[windowLength] = 0.0f;
    }
    else {
      windowBuffer[0] = 1.0f;
      for (i=1; i<=windowLength; i++)
        windowBuffer[i] = (float)(0.54 + 0.46*cos(ftmp*(double)i));
    }
}

void daicsp::Vonhann(float* windowBuffer, int windowLength, int even) 
{
    float ftmp;
    int i;

    ftmp = PI_F/windowLength;

    if (even) {
      for (i=0; i<windowLength; i++)
        windowBuffer[i] = (float)(0.5 + 0.5 * cos(ftmp*((double)i+0.5)));
      windowBuffer[windowLength] = 0.0f;
    }
    else {
      windowBuffer[0] = 1.0f;
      for (i=1; i<=windowLength; i++)
        windowBuffer[i] = (float)(0.5 + 0.5 * cos(ftmp*(double)i));
    }
}

double daicsp::Besseli(double x)
{
    double ax, ans;
    double y;

    if (( ax = fabs( x)) < 3.75)     {
      y = x / 3.75;
      y *= y;
      ans = (1.0 + y * ( 3.5156229 +
                         y * ( 3.0899424 +
                               y * ( 1.2067492 +
                                     y * ( 0.2659732 +
                                           y * ( 0.360768e-1 +
                                                 y * 0.45813e-2))))));
    }
    else {
      y = 3.75 / ax;
      ans = ((exp ( ax) / sqrt(ax))
             * (0.39894228 +
                y * (0.1328592e-1 +
                     y * (0.225319e-2 +
                          y * (-0.157565e-2 +
                               y * (0.916281e-2 +
                                    y * (-0.2057706e-1 +
                                         y * (0.2635537e-1 +
                                              y * (-0.1647633e-1 +
                                                   y * 0.392377e-2)))))))));
    }
    return ans;
}