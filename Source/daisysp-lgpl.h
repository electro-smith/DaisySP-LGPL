#pragma once
#ifndef DSYSP_LGPL_H
#define DSYSP_LGPL_H

/** Control Modules */
#include "Control/line.h"

/** Dynamics Modules */
#include "Dynamics/balance.h"
#include "Dynamics/compressor.h"

/** Effects Modules */
#include "Effects/bitcrush.h"
#include "Effects/fold.h"
#include "Effects/reverbsc.h"

/** Filter Modules */
#include "Filters/allpass.h"
#include "Filters/atone.h"
#include "Filters/biquad.h"
#include "Filters/comb.h"
#include "Filters/mode.h"
#include "Filters/moogladder.h"
#include "Filters/nlfilt.h"
#include "Filters/tone.h"

/** Physical Modeling Modules */
#include "PhysicalModeling/pluck.h"
#include "PhysicalModeling/PolyPluck.h"

/** Synthesis Modules */
#include "Synthesis/blosc.h"

/** Utility Modules */
#include "Utility/jitter.h"
#include "Utility/port.h"

#endif