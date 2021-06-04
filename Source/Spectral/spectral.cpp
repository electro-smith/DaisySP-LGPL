#include <math.h>
#include "dsp.h"
#include "spectral.h"

using namespace daicsp;

float mod2Pi(float value) {
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