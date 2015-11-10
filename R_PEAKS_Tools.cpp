#include "R_PEAKS_Tools.h"

#include <algorithm>

void Tools::convolve(const std::vector<float>& signal, const std::vector<float>& kernel, std::vector<float>& result) const
{
	size_t signalLen = signal.size();
	size_t kernelLen = kernel.size();
	size_t kmin = 0;
	size_t kmax = 0;
	for (size_t i = 0; i < signalLen + kernelLen - 1; ++i)
	{
		result[i] = 0;

		kmin = (i >= kernelLen - 1) ? i - (kernelLen - 1) : 0;
		kmax = (i < signalLen - 1) ? i : signalLen - 1;

		for (size_t j = kmin; j <= kmax; ++j)
		{
			result[i] += signal[j] * kernel[i - j];
		}
	}
}

void Tools::vector_square(std::vector<float>& vec) const
{
	for (size_t i = 0; i < vec.size(); ++i)
	{
		vec[i] = pow(vec[i], 2);
	}
}

void Tools::adjust(std::vector<unsigned int>& peaks, const std::vector<float>& ecgData, float samplingFrequency) const
{
	size_t m = ceil(0.05*samplingFrequency);
	size_t i = 0;

	if(peaks.size() == 0) return;

	if(peaks[i] < m)
	{
		if(peaks[i] + m > ecgData.size())
		{
			auto pos = std::max_element(ecgData.begin(), ecgData.end());
			peaks[i] = pos - ecgData.begin();
		}
		else
		{
			auto pos = std::max_element(ecgData.begin(), ecgData.begin() + peaks[i]);
			peaks[i] = pos - ecgData.begin();
		}

	}
	else
	{
		if(peaks[i] + m > ecgData.size())
		{
			auto pos = std::max_element(ecgData.begin() + (peaks[i] - m), ecgData.end());
			peaks[i] = peaks[i] - m + (pos - (ecgData.begin() + (peaks[i] - m)));
		}
		else
		{
			auto pos = std::max_element(ecgData.begin() + (peaks[i] - m), ecgData.begin() + peaks[i] + m);
			peaks[i] = peaks[i] - m + (pos - (ecgData.begin() + (peaks[i] - m)));
		}
	}

	for(i = 1; i < peaks.size() - 1; ++i)
	{
		auto pos = std::max_element(ecgData.begin() + (peaks[i] - m), ecgData.begin() + peaks[i] + m);
		peaks[i] = peaks[i] - m + (pos - (ecgData.begin() + (peaks[i] - m)));
	}

	i = peaks.size() - 1;
	if(peaks[i] + m > ecgData.size())
	{
		auto pos = std::max_element(ecgData.begin() + (peaks[i] - m), ecgData.end());
		peaks[i] = peaks[i] - m + (pos - (ecgData.begin() + (peaks[i] - m)));
	}
	else
	{
		auto pos = std::max_element(ecgData.begin() + (peaks[i] - m), ecgData.begin() + peaks[i] + m);
		peaks[i] = peaks[i] - m + (pos - (ecgData.begin() + (peaks[i] - m)));
	}
}


void Tools::createTimeVec(std::vector<float>& time, float samplingFrequency, size_t size) const
{
	time.reserve(size);
	for(float i = 0; i < size - 1; ++i)
	{
		time.push_back(i/samplingFrequency);
	}
}
