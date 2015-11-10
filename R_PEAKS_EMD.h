#pragma once

#include "R_PEAKS_IAlgorithm.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>

typedef std::vector<double>  d_vect;

class EMD : public Algorithm
{
  public:
    void process(std::vector<float> const & ecgData, std::vector<unsigned int>& output, float samplingFrequency) const override;
  private:
    void process(const gsl_vector* x, gsl_vector* imf) const;
    void diff(gsl_vector* dx, const gsl_vector* x) const;
    double sum(const gsl_vector* v) const;
    void divide_maxmin(const gsl_vector* analyzedSignal, const d_vect& maxmin, d_vect& mins, d_vect& mins_val, d_vect& maxes, d_vect& maxes_val) const;
    void find_maxmin(gsl_vector* diffs, d_vect& maxmin) const;
    int sign(double val) const;
};
