#pragma once

#include "R_PEAKS_IAlgorithm.h"

class PanTompkins: public Algorithm
{
public:
    void process(std::vector<float> const & ecgData,
                 std::vector<unsigned int>& output,
                 float samplingFrequency) const override;
private:
    void threshold(std::vector<float> const & ecgData,
                   std::vector<float>& signal,
                   float samplingFrequency,
                   std::vector<unsigned int>& output) const;
};
