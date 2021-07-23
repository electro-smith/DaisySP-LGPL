//    DaiCsP is a DSP Library targeted at the Electrosmith Daisy Product Line.
//    Author: Stephen Hensley, 2019
//
//    However, this is decoupled from the hardware in such a way that it
//        should be useful outside of the ARM context with different build configurations.
//
//    A few general notes about the contents of the library:
//        - all memory usage is static.
//        - in cases of potentially large memory usage: the user will either supply a buffer and a size, or the class will be a template that can have size set at compile time.
//        - all modules will have an Init() function, and a Process() function.
//        - all modules, unless otherwise noted, will process a single sample at a time.
//        - all processing will be done with 'float' type unless otherwise noted.
//

#pragma once
#ifndef DCSP_H
#define DCSP_H

/** Control Modules */
#include "Control/line.h"
#include "Control/phasor.h"

/** Dynamics Modules */
#include "Dynamics/balance.h"

/** Effects Modules */
#include "Effects/fold.h"
#include "Effects/reverbsc.h"

/** Filter Modules */
#include "Filters/allpass.h"
#include "Filters/atone.h"
#include "Filters/biquad.h"
#include "Filters/comb.h"
// #include "Filters/mode.h"
#include "Filters/moogladder.h"
#include "Filters/nlfilt.h"
#include "Filters/tone.h"

/** Physical Modeling Modules */
#include "PhysicalModeling/pluck.h"

/** Spectral Modules */
#include "Spectral/spectral.h"
#include "Spectral/spectralanalyzer.h"
#include "Spectral/phasevocoder.h"
#include "Spectral/spectralblur.h"
#include "Spectral/spectralfreeze.h"
#include "Spectral/spectralscale.h"
#include "Spectral/spectralsmooth.h"

/** Utility Modules */
#include "Utility/dsy_lgpl_dsp.h"
#include "Utility/jitter.h"
#include "Utility/port.h"

#endif
