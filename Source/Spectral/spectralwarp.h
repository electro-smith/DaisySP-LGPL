#pragma once
#ifndef DSY_SPECTRALWARP_H
#define DSY_SPECTRALWARP_H

namespace daicsp
{
/** SpectralWarp
 *  Author: Gabriel Ball
 *  Date: 2021-06-14
 * 
 *  ported from Csound pvswarp
 */
class SpectralWarp
{
  public:
    SpectralWarp() {}
    ~SpectralWarp() {}

    /** Initializes the SpectralWarp module.
         *  \param p - description
         */
    void Init();

    /** Processes a single sample and returns it.
         *  \param in - input sample
         */
    float Process(const float &in);

    /** Setter for private member param
         */
    inline void SetParam(const float &param) { param_ = param; }

    /** Setter for private member complex_param
         */
    void SetComplexParam(const float &complex_param);

  private:
    float param_;
    float foo_bar_;
    float a_, b_;
};

} // namespace daicsp

#endif // DSY_SPECTRALWARP_H
