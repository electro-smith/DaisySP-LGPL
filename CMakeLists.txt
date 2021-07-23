
cmake_minimum_required(VERSION 3.15)

project(DAICSP VERSION 0.0.1)

set(SOURCE_FILES 
Source/Control/line.cpp
Source/Control/phasor.cpp
Source/Dynamics/balance.cpp
Source/Effects/fold.cpp
Source/Effects/reverbsc.cpp
Source/Filters/allpass.cpp
Source/Filters/atone.cpp
Source/Filters/biquad.cpp
Source/Filters/comb.cpp
Source/Filters/mode.cpp
Source/Filters/moogladder.cpp
Source/Filters/nlfilt.cpp
Source/Filters/tone.cpp
Source/PhysicalModeling/pluck.cpp
Source/Utility/jitter.cpp
Source/Utility/port.cpp
Source/Spectral/spectral.cpp
Source/Spectral/spectralanalyzer.cpp
Source/Spectral/phasevocoder.cpp
Source/Spectral/spectralblur.cpp
Source/Spectral/spectralfreeze.cpp
Source/Spectral/spectralscale.cpp
Source/Spectral/spectralsmooth.cpp
Source/Spectral/spectralshift.cpp
)

set( HEADER_FILES
Source/Control/line.h
Source/Control/phasor.h
Source/Dynamics/balance.h
Source/Effects/fold.h
Source/Effects/reverbsc.h
Source/Filters/allpass.h
Source/Filters/atone.h
Source/Filters/biquad.h
Source/Filters/comb.h
Source/Filters/mode.h
Source/Filters/moogladder.h
Source/Filters/nlfilt.h
Source/Filters/tone.h
Source/PhysicalModeling/pluck.h
Source/Utility/dsp.h
Source/Utility/jitter.h
Source/Utility/port.h
Source/Spectral/spectral.h
Source/Spectral/spectralanalyzer.h
Source/Spectral/phasevocoder.h
Source/Spectral/spectralblur.h
Source/Spectral/spectralfreeze.h
Source/Spectral/spectralscale.h
Source/Spectral/spectralsmooth.h
Source/Spectral/spectralshift.h
)

include_directories(Source)
include_directories(Source/Control)
include_directories(Source/Dynamics)
include_directories(Source/Effects)
include_directories(Source/Filters)
include_directories(Source/PhysicalModeling)
include_directories(Source/Spectral)
include_directories(Source/Utility)

set(CMAKE_CXX_STANDARD 14)

ADD_LIBRARY( DaiCsP STATIC 
    ${SOURCE_FILES} 
    ${HEADER_FILES} 
    Source/daicsp.h )
    TARGET_INCLUDE_DIRECTORIES(DaiCsP PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Source>" "$<INSTALL_INTERFACE:Source>")
