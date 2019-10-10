#include <algorithm>
#include "binary_segmentation.h"
#pragma once

class Contrast {
public:
    int tL; // left endpt in contrast vector
    int tR; // right endpt in contrast vector
    int thj; // center of contrast
    double nuTnu_inv; // 1 / sum(nu ^ 2)
    double nuTy;
    double shift;
    Contrast(int l, int t, int r, const VecDouble & Sy);
    Contrast();
    void update_nuTY(const VecDouble & y);
    void set_shift(double a);
};

class IntervalInspection {
public:
    int start;
    int end;
    int tau;
    int change_sign;
    IntervalInspection(int s, int e, int t, int d);
};

class PhiInterval {
public:
    double L;
    double U;
    int contained;
    PhiInterval(double L, double U, int contained);
    PhiInterval();
    void print();
};

class BSInferenceOpts {
public:
    int type;
    int n_steps;
    double lower_threshold;
    double upper_threshold;
    double delta;
    int maxit;
    double TOL;
    BSInferenceOpts(int t, int n, double l, double u, double d, int m, double TOL);
    BSInferenceOpts();
};

BSInferenceOpts BSInferenceOpts_from_simple(int t, int n);

class BSInference: public BinarySegmentation {
public:
    BSInferenceOpts opts;
    VecInt sorted_changepoints;
    VecDouble pvals;
    VecDouble approximation_errors;
    std::vector<std::vector<PhiInterval>> phi_intervals;
    void print();
    void add_changepoint_inference(double pval, double approximation_error, std::vector<PhiInterval> phi_interval); // same order as sorted changepoints!
    BSInference(VecInt ordered_changepoints, VecInt changepoint_signs, BSInferenceOpts o, int n):BinarySegmentation(ordered_changepoints, changepoint_signs) {
      opts = o;
      pvals = VecDouble();
      phi_intervals = std::vector<std::vector<PhiInterval>>();
      VecInt s_changepoints(ordered_changepoints);
      s_changepoints.emplace_back(-1);
      s_changepoints.emplace_back(n - 1);
      std::sort(s_changepoints.begin(), s_changepoints.end());
      sorted_changepoints = s_changepoints;
    };
};

BSInference k_step_bs_inference(const VecDouble & y, const int n_steps, const int type, const double sigma, const int window_size, const double approximation_threshold);

std::vector<PhiInterval> partition_phi_space(const VecDouble & y, Contrast & nu, const BSInferenceOpts & opts,
                                             const BinarySegmentation & orig_segmentation);
bool partition_sort(PhiInterval a, PhiInterval b);