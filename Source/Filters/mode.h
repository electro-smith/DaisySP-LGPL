/*
Copyright (c) 2023 Electrosmith, Corp, Francois Blanc, Steven Yi, Paul Batchelor

Use of this source code is governed by the LGPL V2.1
license that can be found in the LICENSE file or at
https://opensource.org/license/lgpl-2-1/
*/

#pragma once
#ifndef DAISY_MODE
#define DAISY_MODE

namespace daisysp
{
/** Resonant Modal Filter */
class Mode
{
  public:
    Mode() {}
    ~Mode() {}
    /** Initializes the instance of the module.
        sample_rate: frequency of the audio engine in Hz
    */
    void Init(float sample_rate);

    /** Processes one input sample through the filter, and returns the output.
    */
    float Process(float in);

    /** Clears the filter, returning the output to 0.0
    */
    void Clear();

    /** Sets the resonant frequency of the modal filter.
        Range: Any frequency such that sample_rate / freq < PI (about 15.2kHz at 48kHz)
    */
    inline void SetFreq(float freq) { freq_ = freq; }
    /** Sets the quality factor of the filter.
        Range: Positive Numbers (Good values range from 70 to 1400)
    */
    inline void SetQ(float q) { q_ = q; }

  private:
    float freq_, q_;
    float xnm1_, ynm1_, ynm2_, a0_, a1_, a2_;
    float d_, lfq_, lq_, sr_;
};
} // namespace daisysp

#endif
