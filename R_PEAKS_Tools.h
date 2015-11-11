#pragma once
#include <vector>
#include <stddef.h>

class Tools
{
  public:
    void convolve(const std::vector<float>& signal, const std::vector<float>& kernel, std::vector<float>& result) const;
    void vector_square(std::vector<float>& vec) const;
    void adjust(std::vector<unsigned int>& peaks, const std::vector<float>& ecg, float samplingFrequency) const;
    void createTimeVec(std::vector<float>& time, float samplingFrequency, size_t size) const;
    void writeVectorToFile(std::vector<float> signal, char * path) const;
    void writeVectorToFile(std::vector<unsigned int> signal, char * path) const;
};


