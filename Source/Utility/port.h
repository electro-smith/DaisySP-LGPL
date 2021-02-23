
#pragma once
#ifndef DSY_PORT_H
#define DSY_PORT_H
#ifdef __cplusplus

namespace daicsp
{
/** Applies portamento to an input signal. 

At each new step value, the input is low-pass filtered to 
move towards that value at a rate determined by ihtim. ihtim is the half-time of the 
function (in seconds), during which the curve will traverse half the distance towards the new value, 
then half as much again, etc., theoretically never reaching its asymptote.

Original Author(s): Robbin Whittle, John ffitch

Year: 1995, 1998

Location: Opcodes/biquad.c

Ported by: Stephen Hensley
*/
class Port
{
  public:
    Port() {}
    ~Port() {}
    /** Initializes Port module

        \param sample_rate: sample rate of audio engine
        \param htime: half-time of the function, in seconds.
    */

    void Init(float sample_rate, float htime);

    /** Applies portamento to input signal and returns processed signal. 
        \return slewed output signal
    */
    float Process(float in);


    /** Sets htime
    */
    inline void SetHtime(float htime) { htime_ = htime; }
    /** returns current value of htime
    */
    inline float GetHtime() { return htime_; }

  private:
    float htime_;
    float c1_, c2_, yt1_, prvhtim_;
    float sample_rate_, onedsr_;
};
} // namespace daicsp
#endif
#endif
