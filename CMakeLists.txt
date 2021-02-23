
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
)

include_directories(${CMAKE_SOURCE_DIR}/Source)
include_directories(${CMAKE_SOURCE_DIR}/Source/Control)
include_directories(${CMAKE_SOURCE_DIR}/Source/Dynamics)
include_directories(${CMAKE_SOURCE_DIR}/Source/Effects)
include_directories(${CMAKE_SOURCE_DIR}/Source/Filters)
include_directories(${CMAKE_SOURCE_DIR}/Source/PhysicalModeling)
include_directories(${CMAKE_SOURCE_DIR}/Source/Utility)

set(CMAKE_CXX_STANDARD 14)

ADD_LIBRARY( DaisySP STATIC 
    ${SOURCE_FILES} 
    ${HEADER_FILES} 
    Source/daisysp.h )

