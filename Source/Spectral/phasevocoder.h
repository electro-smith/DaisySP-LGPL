#pragma once
#ifndef DSY_PHASEVOCODER_H
#define DSY_PHASEVOCODER_H

namespace daicsp
{

/** PhaseVocoder
 *  Author: Gabriel Ball
 *  Date: 2021-06-03
 */
class PhaseVocoder
{
    public:
        PhaseVocoder () {}
        ~PhaseVocoder () {}

        /** Initializes the PhaseVocoder module.
         *  \param p - description
         */
        void Init();

        /** Processes a single sample and returns it.
         *  \param in - input sample
         */
        float Process(const float &in);

    private:
        
};

} // namespace daicsp

#endif // DSY_PHASEVOCODER_H
