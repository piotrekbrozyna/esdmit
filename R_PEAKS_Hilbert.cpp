#include "R_PEAKS_Hilbert.h"

#include <cmath>
#include <algorithm>
#include <inttypes.h>
#include <fftw3.h>

#include <iostream>


void Hilbert::process(std::vector<float> const & ecgData,
        std::vector<unsigned int>& output, float samplingFrequency) const
{
    std::vector<float> kernel = { -0.125, -0.250, 0, 0.250, 0.125 };
    std::vector<float> y(ecgData.size() + kernel.size() - 1);
    std::vector<float> y2;

    m_tools.writeVectorToFile(ecgData, "HilbertInput.csv", true);
    
    // differentiation
    m_tools.convolve(ecgData, kernel, y);
    y = std::vector<float>(y.begin()+2,y.end()-2);

    y2.reserve(y.size());
    hilbert(y,y2);

    for (size_t i = 0; i < y2.size(); ++i)
        y2[i] = sqrt(y[i]*y[i] + y2[i]*y2[i]);

    size_t k = 5 * 60 * samplingFrequency;
    size_t size = ecgData.size();

    std::vector<unsigned int> o;
    o.reserve(k);
    for (size_t i = 0; i < size; i += k)
    {
        o.clear();
        if (i + k >= size)
        {
            threshold(
                    std::vector<float> (ecgData.begin() + i, ecgData.end()),
                    std::vector<float> (y2.begin() + i, y2.end()),
                    samplingFrequency, o);
        }
        else
        {
            threshold(
                    std::vector<float> (ecgData.begin() + i, ecgData.begin() + i + k),
                    std::vector<float> (y2.begin() + i, y2.begin() + i + k),
                    samplingFrequency, o);
        }

        for(size_t j = 0; j < o.size(); ++j)
        {
            output.push_back(o[j] + i);
        }
    }

    m_tools.writeVectorToFile(output, "HilbertOutput.csv", true);
}

void Hilbert::process_on_imf(std::vector<float> const & ecgData, std::vector<float>& imf,
        std::vector<unsigned int>& output, float samplingFrequency) const
{
    std::vector<float> kernel = { -0.125, -0.250, 0, 0.250, 0.125 };
    std::vector<float> y(imf.size() + kernel.size() - 1);
    std::vector<float> y2;

    m_tools.convolve(imf, kernel, y);
    y = std::vector<float>(y.begin()+2,y.end()-2);

    y2.reserve(y.size());
    hilbert(y,y2);

    for (size_t i = 0; i < y2.size(); ++i)
        y2[i] = sqrt(y[i]*y[i] + y2[i]*y2[i]);

    size_t k = 2 * 60 * samplingFrequency;
    size_t size = ecgData.size();

    std::vector<unsigned int> o;
    o.reserve(k);
    for (size_t i = 0; i < size; i += k)
    {
        o.clear();
        if (i + k >= size)
        {
            threshold(
                    std::vector<float> (ecgData.begin() + i, ecgData.end()),
                    std::vector<float> (y2.begin() + i, y2.end()),
                    samplingFrequency, o);
        }
        else
        {
            threshold(
                    std::vector<float> (ecgData.begin() + i, ecgData.begin() + i + k),
                    std::vector<float> (y2.begin() + i, y2.begin() + i + k),
                    samplingFrequency, o);
        }

        for(size_t j = 0; j < o.size(); ++j)
        {
            output.push_back(o[j] + i);
        }
    }
}

void Hilbert::threshold(std::vector<float> ecgData, std::vector<float> y, float samplingFrequency, std::vector<unsigned int>& output) const
{
    float max = *(std::max_element(y.begin(), y.begin()+5*samplingFrequency));

    float threshold = 0.8 * max;

    output.push_back(0);

    float v = std::floor(0.7*samplingFrequency);
    std::vector<float> RR(8,v);
    std::vector<float> tmpECG(5,0.0);

    float mR = 0.0;

    for (size_t i = 0; i < ecgData.size(); ++i)
    {
        float dt = (i-output.back());
        mR = std::accumulate(RR.begin(), RR.end(), 0.0)/8.0;

        if(ecgData[i] > threshold)
        {
            if(dt <= 0.2*samplingFrequency or dt <= 0.55*mR)
            {
                if(ecgData[i] > ecgData[output.back()])
                {
                    output.back() = i;
                    tmpECG[output.size() % 5] = ecgData[output.back()];
                }

                if (output.size() > 1)
                {
                	RR[output.size() % 8] = output.back() - output[output.size()-2];
                }
            }
            else
            {
                output.push_back(i);
                tmpECG[output.size() % 5] = ecgData[output.back()];
                if(output.size() > 1)
                {
                	RR[output.size() % 8] = (output.back() - output[output.size()-2]);
                }
            }
        }
        if(output.size() > 4)
        {
            threshold = 0.55*(std::accumulate(tmpECG.begin(), tmpECG.end(), 0.0)
                            - *(std::max_element(tmpECG.begin(), tmpECG.end()))) / 4.0;
        }
    }

    if(output[0] == 0)
    	output.erase(output.begin());

    m_tools.adjust(output, ecgData, samplingFrequency);
}

void Hilbert::hilbert(std::vector<float>& ecgData, std::vector<float>& output) const
{
    size_t resultSize = ecgData.size();

    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * resultSize);
    fftw_complex *in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * resultSize);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * resultSize);

    for(size_t i = 0; i < resultSize; i++){
       in[i][0] = ecgData[i];
       in[i][1] = 0;
    }

    //FFT
    fftw_plan pf = fftw_plan_dft_1d(resultSize, in, in2, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(pf);
    fftw_destroy_plan(pf);

    for (int i = 0; i < resultSize; i++)
    {
        if (i == 0 or (i == resultSize / 2))
        {
            //in2[i][0] = 0;
            //in2[i][1] = 0;
        }
        else if (i < resultSize / 2)
        {
            in2[i][0] = 2 * in2[i][0];
            in2[i][1] = 2 * in2[i][1];
        }
        else if (i > resultSize / 2)
        {
            in2[i][0] = 0 * in2[i][0];
            in2[i][1] = 0 * in2[i][1];
        }
    }

    //IFFT
    fftw_plan pb = fftw_plan_dft_1d(resultSize, in2, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pb);
    fftw_destroy_plan(pb);

    output.reserve(resultSize);
    for(size_t i = 0; i < resultSize; i++){
       output.push_back(out[i][0]/resultSize);
    }

    fftw_free(in);
    fftw_free(in2);
    fftw_free(out);
}
