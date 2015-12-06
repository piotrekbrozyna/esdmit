#include "R_PEAKS_PanTompkins.h"

#include <cmath>
#include <inttypes.h>
#include <algorithm>
#include <iostream>

void PanTompkins::process(std::vector<float> const & ecgData, std::vector<unsigned int>& output, float samplingFrequency) const
{
    std::vector<float> kernel = { -0.125, -0.250, 0, 0.250, 0.125 };

    std::vector<float> y(ecgData.size() + kernel.size() - 1);
    
    m_tools.writeVectorToFile(ecgData, "Input.csv", true);

    // derivation
    m_tools.convolve(ecgData, kernel, y);

    m_tools.vector_square(y);

    float n = float(150) / 1000;

    uint32_t N = floor(samplingFrequency * n);

    std::vector<float> kernel2(N, float(1) / N);

    std::vector<float> y1(y.size() + kernel2.size() - 1);

    // integrating
    m_tools.convolve(y, kernel2, y1);
    //m_tools.writeVectorToFile(y1, "AdditionalData.csv", true);

    size_t dl = y1.size() - ecgData.size();

    std::vector<float> y3(y1.begin() + (dl / 2), y1.end() - (dl / 2) - dl % 2);

    threshold(ecgData, y3, samplingFrequency, output);
    m_tools.writeVectorToFile(output, "PanTompkinsOutput.csv", true);
    //m_tools.adjust(output, ecgData, samplingFrequency);
}

void PanTompkins::threshold(std::vector<float> const & ecgData, std::vector<float>& signal, float samplingFrequency, std::vector<unsigned int>& output) const
{
<<<<<<< HEAD
    std::vector<unsigned int> fidualMark;
    fidualMark.reserve(ecgData.size());

    std::vector<float> kernel = { -1.0, 1.0 };
    std::vector<float> dy(signal.size() + kernel.size() - 1);

    for(auto&& x : signal)
        x = (x <= 0.0001) ? 0.0 : x;

    // rising edge detection
    m_tools.convolve(signal, kernel, dy);

    for(unsigned int i = 2; i < dy.size()-1; ++i)
    {
        if((dy[i-1] > 0) && (dy[i] <= 0))
        {
            if(fidualMark.size() > 0)
            {
                // minimal distance between consecutive points - 200 ms
                if(i - fidualMark.back() >= ceil(0.2*samplingFrequency))
                {
                    fidualMark.push_back(i);
                }
                else
                {
                    if(signal[i] >= signal[fidualMark.back()])
                    {
                        fidualMark.back() = i;
                    }
                }
            }
            else
                fidualMark.push_back(i);
        }
    }

    /*
    filtered signal:
    - 
    */
    // signal after filtration - spkf, npkf, thf1, thf2
    // signal after integrating - spki, npki, thi1, thi2

    // maximum from first 2 seconds of integrated signal??
    float spki = *(std::max_element(signal.begin(), signal.begin() + 2*samplingFrequency))/3.0;
    float npki = std::accumulate(signal.begin(), signal.begin() + 2*samplingFrequency,0.0)/(2*samplingFrequency*2.0);
    float thi1 = spki;
    float thi2 = npki;
    float spkf = *(std::max_element(ecgData.begin(), ecgData.begin() + 2*samplingFrequency))/3.0;
    float npkf = std::accumulate(ecgData.begin(), ecgData.begin() + 2*samplingFrequency,0.0)/(2*samplingFrequency*2.0);
    float thf1 = spkf;
    float thf2 = npkf;

    // window size (neighborhood) for finding local maxima
    size_t n = ceil(0.1*samplingFrequency);
    // if true, peak is a T wave and should be skipped
    bool skip = false;
    // avarage R-to-R interval
    float mRR = 0;
    size_t pos = 0;
    size_t posTmp = 0;
    // R-to-R inetrvals
    std::vector<unsigned int> RR;
    RR.reserve(8);
    unsigned int fidualMarkTmp = 0;
    // slopes of 2 consecutive peaks
    float slope1, slope2;

    for(size_t i = 0; i < fidualMark.size(); ++i)
    {
        if((fidualMark[i] > n) && (fidualMark[i] + n < ecgData.size()))
        {
            pos = std::max_element(ecgData.begin() + fidualMark[i] - n, ecgData.begin() + fidualMark[i] + n) - ecgData.begin();
        }
        else if(i == 0)
        {
            pos = std::max_element(ecgData.begin(), ecgData.begin() + fidualMark[i] + n) - ecgData.begin();
        }
        else if(fidualMark[i] + n > ecgData.size())
        {
            pos = std::max_element(ecgData.begin() + fidualMark[i] - n, ecgData.end()) - ecgData.begin();
        }

        if (output.size() > 8)
        {
            // compute last 8 R-to-R intervals
            for(int i = 0; i < 8; ++i)
                RR[8-i] = *(output.end() - i - 1) - *(output.end() - i - 2);

            // compute average R-to-R interval
            mRR = std::accumulate(RR.begin(), RR.end(), 0.0) / RR.size();

            // uneven heart pulse? then reduce thresholds by half
            if((RR.back() <= 0.92*mRR) || (RR.back() >= 1.16*mRR))
            {
                thf1 = 0.5*thf1;
                thf2 = 0.5*thf1;
            }
        }

        if(mRR > 0)
        {
            // R peak omitted? then search back
            if(fidualMark[i] - output.back() >= 1.66*mRR)
            {
                fidualMarkTmp = std::max_element(signal.begin() + output.back() + floor(0.2*samplingFrequency), signal.begin() + fidualMark[i] - floor(0.2*samplingFrequency)) - (signal.begin() + output.back() + floor(0.2*samplingFrequency));
                // WTF? adding what was subtracted in previous step?
                fidualMarkTmp = fidualMarkTmp + output.back() + floor(0.2*samplingFrequency);

                if((fidualMarkTmp > n) && (fidualMarkTmp + n < ecgData.size()))
                {
                    posTmp = std::max_element(ecgData.begin() + fidualMarkTmp - n, ecgData.begin() + fidualMarkTmp + n) - (ecgData.begin() + fidualMarkTmp - n);
                    // WTF no. 2
                    posTmp = fidualMarkTmp - n + posTmp;
                }
                else if(i == 0)
                {
                    posTmp = std::max_element(ecgData.begin(), ecgData.begin() + fidualMarkTmp + n) - ecgData.begin();
                }
                else if(fidualMarkTmp + n > ecgData.size())
                {
                    posTmp = std::max_element(ecgData.begin() + fidualMarkTmp - n, ecgData.end()) - (ecgData.begin() + fidualMarkTmp - n);
                    // WTF no. 3
                    posTmp = fidualMarkTmp - n  + posTmp;
                }

                if(signal[fidualMarkTmp] > thi2)
                {
                    if(ecgData[posTmp] > thf2)
                    {
                        output.push_back(posTmp);
                        spkf = 0.25*ecgData[posTmp] + 0.75*spkf;
                    }

                    spki = 0.25*signal[fidualMarkTmp] + 0.75*spki;
                }
            }
        }

        if(signal[fidualMark[i]] >= thi1)
        {
            if(output.size() > 2)
            {
                // peak can be a T wave
                if(fidualMark[i] - output.back() <= ceil(0.360*samplingFrequency))
                {
                    slope1 = signal[fidualMark[i]-4]*-1 + signal[fidualMark[i]-3]*-2 + signal[fidualMark[i]-1]*2 + signal[fidualMark[i]];
                    slope2 = signal[output.back()-4]*-1 + signal[output.back()-3]*-2 + signal[output.back()-1]*2 + signal[output.back()];

                    if(abs(slope1) <= 0.5*abs(slope2))
                    {
                        skip = true;
                        npki = 0.125*signal[fidualMark[i]] + 0.875*npki;
                        npkf = 0.125*ecgData[pos] + 0.875*npkf;
                    }
                    else
                        skip = false;
                }
            }

            if(!skip)
            {
                if(ecgData[pos] >= thf1)
                {
                    output.push_back(pos);
                    spkf = 0.125*ecgData[pos] + 0.875*spkf;
                }

                spki = 0.125*signal[fidualMark[i]] + 0.875*spki;
            }
        }
        else if((thi2 <= signal[fidualMark[i]] && signal[fidualMark[i]] < thi1) || signal[fidualMark[i]] < thi2)
        {
            npki = 0.125*signal[fidualMark[i]] + 0.875*npki;
            npkf = 0.125*ecgData[pos] + 0.875*npkf;
        }


        thi1 = npki + 0.25*(spki-npki);
        thi2 = 0.5*thi1;
        thf1 = npkf + 0.25*(spkf-npkf);
        thf2 = 0.5*thf1;

        skip = false;
    }
    m_tools.writeVectorToFile(fidualMark, "AdditionalData.csv", false);
    std::vector<float> finalThreshold = { thi1, thi2, thf1, thf2 };
    std::vector<float> cos= { spki, spkf, npki, npkf };
    m_tools.writeVectorToFile(finalThreshold, "AdditionalData.csv", false);
    m_tools.writeVectorToFile(cos, "AdditionalData.csv", false);
}
