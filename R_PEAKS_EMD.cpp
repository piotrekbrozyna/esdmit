#include "R_PEAKS_EMD.h"
#include "R_PEAKS_Hilbert.h"
#include "R_PEAKS_PanTompkins.h"

#include <iostream>
#include <cmath>


void EMD::process(std::vector<float> const & ecgData, std::vector<unsigned int>& output, float samplingFrequency) const
{
    //m_tools.writeVectorToFile(ecgData, "EMDInput.csv", true);
    gsl_vector* ekg = gsl_vector_alloc(ecgData.size());
    for (int i = 0; i < ecgData.size(); i++)
    {
        gsl_vector_set(ekg, i, (double)ecgData[i]);
    }

    gsl_vector* _imf = gsl_vector_alloc(ecgData.size());

    process(ekg, _imf);

    std::vector<float> imf;
    imf.reserve(_imf->size);
    for (int i = 0; i < _imf->size; i++)
    {
        imf.push_back((float)gsl_vector_get(_imf, i));
    }

    Hilbert().process_on_imf(ecgData, imf, output, samplingFrequency);
    m_tools.writeVectorToFile(output, "EMDOutput.csv", true);
}

void EMD::process(const gsl_vector* x, gsl_vector* imf) const
{
    unsigned int length = x->size;
    unsigned int i, j, k, p, N;
    double SD; //standard deviation
    double eps = 0.0000001; //to avoid zero values

    std::vector<double> maxmin; // to store the optima (min and max without distinction so far)
    std::vector<double> maxes, mins, maxes_val, mins_val;

    //spline variables
    gsl_interp_accel* acc_min;
    gsl_spline* spline_min;
    gsl_interp_accel* acc_max;
    gsl_spline* spline_max;

    std::wcout << "Length of EKG signal: " << length << endl;
    int max_parts = ceil(length / 1000.0);
    std::wcout << "Number of parts: " << max_parts << endl;

    unsigned int prog;

    for (p = 0; p < max_parts * 1000; p += 1000)
    {
        if (p == (max_parts - 1) * 1000) // last part of signal
            N = length - p;
        else
            N = 1000;
        if(N==0)
            continue;
        //std::cout << N << endl;
        gsl_vector* x_part = gsl_vector_alloc(N);
        for (i = 0; i < N; ++i)
        {
            gsl_vector_set(x_part, i, gsl_vector_get(x, p + i));
        }

        gsl_vector* c = gsl_vector_alloc(N);
        gsl_vector* h = gsl_vector_alloc(N);
        gsl_vector* prevh = gsl_vector_alloc(N);
        gsl_vector* prevh_temp = gsl_vector_alloc(N);
        gsl_vector* d = gsl_vector_alloc(N - 1); //store diffs
        gsl_vector* mean = gsl_vector_alloc(N);
        gsl_vector* maxenv = gsl_vector_alloc(N);
        gsl_vector* minenv = gsl_vector_alloc(N);

        //copy of the input signal
        gsl_vector_memcpy(c, x_part);

        //inner loop to find each imf
        gsl_vector_memcpy(h, c); // at the beginning of the sifting process, h is the signal
        SD = 1; //Standard deviation which will be used to stop the sifting process

        while (SD > 0.3)
        { // while the standard deviation is higher than 0.3 (typical value)

            //find local max/min points
            maxmin.clear();
            diff(d, h); //approximate derivative
            find_maxmin(d, maxmin);

            if (maxmin.size() < 2)  //then it is the residue
                break;

            //clearing previous values
            maxes.clear();
            mins.clear();
            maxes_val.clear();
            mins_val.clear();

            //divide maxmin into maxes and mins
            divide_maxmin(h, maxmin, mins, mins_val, maxes, maxes_val);

            //make endpoints both maxes and mins
            maxes.insert(maxes.begin(), 0);
            maxes.push_back(N - 1);
            maxes_val.insert(maxes_val.begin(), gsl_vector_get(h, 0));
            maxes_val.push_back(gsl_vector_get(h, N - 1));

            mins.insert(mins.begin(), 0);
            mins.push_back(N - 1);
            mins_val.insert(mins_val.begin(), gsl_vector_get(h, 0));
            mins_val.push_back(gsl_vector_get(h, N - 1));

            //spline interpolate to get max and min envelopes; form imf
            gsl_interp_accel* acc_max = gsl_interp_accel_alloc();
            gsl_interp_accel* acc_min = gsl_interp_accel_alloc();
            spline_min = gsl_spline_alloc(gsl_interp_cspline, mins.size());
            spline_max = gsl_spline_alloc(gsl_interp_cspline, maxes.size());

            gsl_spline_init(spline_min, (const double *) mins.data(),
                    (const double *) mins_val.data(), mins.size());
            gsl_spline_init(spline_max, (const double *) maxes.data(),
                    (const double *) maxes_val.data(), maxes.size());

            const size_t mins_end = mins.at(mins.size() - 1);
            const size_t stride_minenv = minenv->stride;

            for (i = mins[0], k = 0; i <= mins_end; ++i, ++k)
            {
                minenv->data[k * stride_minenv] = gsl_spline_eval(spline_min, i,
                        acc_min);
            }

            const size_t maxes_end = maxes.at(maxes.size() - 1);
            const size_t stride_maxenv = maxenv->stride;

            for (i = maxes[0], k = 0; i <= maxes_end; ++i, ++k)
            {
                maxenv->data[k * stride_maxenv] = gsl_spline_eval(spline_max, i,
                        acc_max);
            }

            gsl_spline_free(spline_min);
            gsl_spline_free(spline_max);
            gsl_interp_accel_free(acc_min);
            gsl_interp_accel_free(acc_max);

            gsl_vector_memcpy(mean, maxenv);
            gsl_vector_add(mean, minenv);
            gsl_vector_scale(mean, 0.5);
            gsl_vector_memcpy(prevh, h); //copy of the previous value of h before modifying it
            gsl_vector_sub(h, mean); //subtract mean to h

            //calculate prevh.^2 + eps
            gsl_vector_memcpy(prevh_temp, prevh);
            gsl_vector_mul(prevh_temp, prevh_temp);
            gsl_vector_add_constant(prevh_temp, eps);

            //calculate ((prevh - h).^2)
            gsl_vector_sub(prevh, h);
            gsl_vector_mul(prevh, prevh);

            gsl_vector_div(prevh, prevh_temp);
            SD = sum(prevh);
        }

        //store the extracted IMF in the matrix imf
        const size_t stride_imf = imf->stride;
        const size_t stride_h = h->stride;

        for (i = 0; i < N; ++i)
        {
            imf->data[(i + p) * stride_imf] = h->data[i * stride_h];
        }

        //stop criterion of the algorithm. if we reach the end before n
        if (maxmin.size() < 2)
        {
            std::cout << "Achieved stop criterion of the algorithm" << std::endl;
            break;
        }

        gsl_vector_sub(c, h);  //subtract the extracted IMF from the signal

        gsl_vector_free(x_part);
        gsl_vector_free(h);
        gsl_vector_free(d);
        gsl_vector_free(c);
        gsl_vector_free(maxenv);
        gsl_vector_free(minenv);
        gsl_vector_free(mean);
        gsl_vector_free(prevh);
        gsl_vector_free(prevh_temp);
    }
}

double EMD::sum(const gsl_vector* v) const
{
    double sum = 0;

    const size_t N = v->size;
    const size_t stride_v = v->stride;

    for (size_t i = 0; i < N; ++i)
    {
        sum += v->data[i * stride_v];
    }

    return sum;
}

void EMD::diff(gsl_vector* dx, const gsl_vector* x) const
{
    const size_t N = dx->size;

    const size_t stride_dx = dx->stride;
    const size_t stride_x = x->stride;

    for (size_t i = 0; i < N; ++i)
    {
        dx->data[i * stride_dx] = x->data[(i + 1) * stride_x] - x->data[i * stride_x];
    }
}

void EMD::find_maxmin(gsl_vector* diffs, d_vect& maxmin) const
{

    const size_t diffs_end = diffs->size - 1;
    const size_t stride_diffs = diffs->stride;

    for (int i = 1; i < diffs_end; ++i) {

        if (diffs->data[i * stride_diffs] == 0) {    //we are on a zero
            if (sign(diffs->data[(i-1) * stride_diffs])
                    != sign(diffs->data[(i+1) * stride_diffs])) { //it is a maximum
                maxmin.push_back(i);
            }
        } else if (sign(diffs->data[i * stride_diffs])
                != sign(diffs->data[(i+1) * stride_diffs])) { //we are straddling a zero so
            maxmin.push_back(i + 1); //define zero as at i+1 (not i)
        }
    }
}

void EMD::divide_maxmin(const gsl_vector* analyzedSignal, const d_vect& maxmin, d_vect& mins,
        d_vect& mins_val, d_vect& maxes, d_vect& maxes_val) const
{
    const size_t stride_sig = analyzedSignal->stride;

    if (maxmin.at(0) > maxmin.at(1))
    { //first one is a max not a min
        for (size_t k = 0; k < maxmin.size(); ++k)
        {
            if (k % 2 == 0)
            {
                maxes.push_back((double) maxmin.at(k));
                maxes_val.push_back(
                        analyzedSignal->data[(size_t) maxmin.at(k) * stride_sig]);
            } else
            {
                mins.push_back((double) maxmin.at(k));
                mins_val.push_back(
                        analyzedSignal->data[(size_t) maxmin.at(k) * stride_sig]);
            }
        }
    }
    else
    {                               //is the other way around
        for (size_t k = 0; k < maxmin.size(); ++k)
        {
            if (k % 2 == 0)
            {
                mins.push_back((double) maxmin.at(k));
                mins_val.push_back(
                        analyzedSignal->data[(size_t) maxmin.at(k) * stride_sig]);
            } else
            {
                maxes.push_back((double) maxmin.at(k));
                maxes_val.push_back(
                        analyzedSignal->data[(size_t) maxmin.at(k) * stride_sig]);
            }
        }
    }
}

inline int EMD::sign(double val) const{
	return (val > 0) - (val < 0);
}
