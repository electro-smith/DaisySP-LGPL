/*
Copyright (c) 2023 Electrosmith, Corp, Victor Lazzarini, John ffitch (fast tanh), Bob Moog

Use of this source code is governed by the LGPL V2.1
license that can be found in the LICENSE file or at
https://opensource.org/license/lgpl-2-1/
*/

#pragma once
#ifndef DSY_MOOGLADDER_H
#define DSY_MOOGLADDER_H

#include <stdint.h>
#ifdef __cplusplus

namespace daisysp
{
/** Moog ladder filter module*/
class MoogLadder
{
  public:
    MoogLadder() {}
    ~MoogLadder() {}
    /** Initializes the MoogLadder module.
        sample_rate - The sample rate of the audio engine being run. 
    */
    void Init(float sample_rate);


    /** Processes the lowpass filter
    */
    float Process(float in);

    /** 
        Sets the cutoff frequency or half-way point of the filter.
        Arguments
        - freq - frequency value in Hz. Range: Any positive value.
    */
    inline void SetFreq(float freq) { freq_ = freq; }
    /** 
        Sets the resonance of the filter.
    */
    inline void SetRes(float res) { res_ = res; }

  private:
    float istor_, res_, freq_, delay_[6], tanhstg_[3], old_freq_, old_res_,
        sample_rate_, old_acr_, old_tune_;
    float my_tanh(float x);
};
} // namespace daisysp
#endif
#endif
