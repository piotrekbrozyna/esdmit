#pragma once

#include <vector>
#include "R_PEAKS_Tools.h"

class Algorithm
{
public:
    virtual void process(std::vector<float> const & ecgData,
                         std::vector<unsigned int>& output,
                         float samplingFrequency) const = 0;

    Algorithm() : m_tools(Tools()) {}
    virtual ~Algorithm() {}

protected:
    Tools m_tools;
};
