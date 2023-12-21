/*
Copyright (c) 2023 Electrosmith, Corp, John FFitch, Gabriel Maldonado

Use of this source code is governed by the LGPL V2.1
license that can be found in the LICENSE file or at
https://opensource.org/license/lgpl-2-1/
*/

#pragma once
#ifndef DSY_FOLD_H
#define DSY_FOLD_H

#include <stdint.h>
#ifdef __cplusplus

namespace daisysp
{
/** Fold module */
class Fold
{
  public:
    Fold() {}
    ~Fold() {}
    /** Initializes the fold module.
    */
    void Init();


    /** applies foldover distortion to input 
    */
    float Process(float in);


    /** 
        \param incr : set fold increment
    */
    inline void SetIncrement(float incr) { incr_ = incr; }

  private:
    float incr_, index_, value_;
    int   sample_index_;
};
} // namespace daisysp
#endif
#endif
