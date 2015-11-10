#pragma once

#include "R_PEAKS_IAlgorithm.h"

class Hilbert : public Algorithm
{
  public:
    void process(std::vector<float> const & ecgData, std::vector<unsigned int>& output, float samplingFrequency) const override;
    void process_on_imf(std::vector<float> const & ecgData, std::vector<float>& imf,
            std::vector<unsigned int>& output, float samplingFrequency) const;
  private:
    void hilbert(std::vector<float>& ecgData, std::vector<float>& output) const;
    void threshold(std::vector<float> ecgData, std::vector<float> y, float samplingFrequency, std::vector<unsigned int>& output) const;
};
