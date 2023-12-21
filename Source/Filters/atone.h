/*
Copyright (c) 2023 Electrosmith, Corp, Barry Vercoe, John FFitch, Gabriel Maldonado

Use of this source code is governed by the LGPL V2.1
license that can be found in the LICENSE file or at
https://opensource.org/license/lgpl-2-1/
*/

#pragma once
#ifndef DSY_ATONE_H
#define DSY_ATONE_H

#include <stdint.h>
#ifdef __cplusplus

namespace daisysp
{
/** A first-order recursive high-pass filter with variable frequency response. */
class ATone
{
  public:
    ATone() {}
    ~ATone() {}
    /** Initializes the ATone module.
        \param sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** Processes one sample through the filter and returns one sample.
        \param in - input signal 
    */
    float Process(float &in);

    /** Sets the cutoff frequency or half-way point of the filter.
        \param freq - frequency value in Hz. Range: Any positive value.
    */
    inline void SetFreq(float &freq)
    {
        freq_ = freq;
        CalculateCoefficients();
    }

    /** get current frequency
        \return the current value for the cutoff frequency or half-way point of the filter.
    */
    inline float GetFreq() { return freq_; }

  private:
    void  CalculateCoefficients();
    float out_, prevout_, in_, freq_, c2_, sample_rate_;
};
} // namespace daisysp
#endif
#endif
