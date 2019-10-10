#include "binary_segmentation.h"
#include <math.h>
#include <numeric>
#include <algorithm>

// class functions

void BinarySegmentation::print(){
  printf("Binary segmentation with %d changepoints\n", static_cast<int>(ordered_changepoints.size()));
  printf("Ordered changepoints:\n");
  for (const auto &i : ordered_changepoints) {
  std::cout << i << '\t';
  }
  std::cout << std::endl;

  printf("Ordered signs:\n");
  for (const auto &i : changepoint_signs) {
    std::cout << i << '\t';
  }
  std::cout << std::endl;
}

BinarySegmentation::BinarySegmentation(const VecInt & o, const VecInt & s) {
  ordered_changepoints = o;
  changepoint_signs = s;
}

BinarySegmentation::BinarySegmentation() {}

CusumInterval::CusumInterval(double m, int a, int s) {
  max_csum = m;
  argmax_csum = a;
  sign_csum = s;
}

CusumInterval::CusumInterval() {
  max_csum = 0;
  argmax_csum = 0;
  sign_csum = 0;
}

// end class defns


double calculate_csum(const VecDouble & Sy, int tau, int s, int e) {
  double right_mean = (Sy[e] - Sy[tau]) / (e - tau);
  double left_mean;
  if (s == 0) {
    left_mean = (Sy[tau] - 0.0) / (tau - s + 1);
  } else {
    left_mean = (Sy[tau] - Sy[s - 1]) / (tau - s + 1);
  }
  double norm_constant = sqrt(1.0 / ( (1.0/ (e - tau)) + (1.0 / (tau + 1 - s)) ));
  return norm_constant * (right_mean - left_mean);
}


CusumInterval argmax_tau_csum_interval(const VecDouble & Sy, int s, int e) {
  double max_csum = -INFINITY;
  double csum = 0;
  int argmax_csum = 0;
  int sign_csum = 0;

  for (int tau = s; tau < e; tau++) {
    csum = calculate_csum(Sy, tau, s, e);
    if (ABS(csum) > max_csum) {
      max_csum = ABS(csum);
      argmax_csum = tau;
      sign_csum = (csum > 0) - (csum < 0);
    }
  }
  return CusumInterval(max_csum, argmax_csum, sign_csum);
}

BinarySegmentation k_step_bs(const VecDouble & Sy, const int n_steps) {
  CusumInterval csums;
  int tau_candidate, sgn_candidate;
  VecInt ordered_changepoints, changepoint_signs, sorted_changepoints;
  sorted_changepoints.emplace_back(-1);
  sorted_changepoints.emplace_back(Sy.size() - 1);
  double seg_max;
  for (int step_i = 0;  step_i < n_steps; step_i++) {
    if (step_i == 0) {
      csums = argmax_tau_csum_interval(Sy, 0, Sy.size() - 1);
      ordered_changepoints.emplace_back(csums.argmax_csum);
      changepoint_signs.emplace_back(csums.sign_csum);
      sorted_changepoints.emplace_back(csums.argmax_csum);
      std::sort(sorted_changepoints.begin(), sorted_changepoints.end());
    } else {
      seg_max = -INFINITY;
      for (int seg_i = 0; seg_i <= step_i; seg_i++) {
        if (sorted_changepoints[seg_i+1] > sorted_changepoints[seg_i] + 1) {
          csums = argmax_tau_csum_interval(Sy, sorted_changepoints[seg_i] + 1, sorted_changepoints[seg_i+1]);
          if (csums.max_csum > seg_max) {
            tau_candidate = csums.argmax_csum;
            sgn_candidate = csums.sign_csum;
            seg_max = csums.max_csum;
          }
        }
      }
      ordered_changepoints.emplace_back(tau_candidate);
      changepoint_signs.emplace_back(sgn_candidate);
      sorted_changepoints.emplace_back(tau_candidate);
      std::sort(sorted_changepoints.begin(), sorted_changepoints.end());
    }
  }
  return BinarySegmentation(ordered_changepoints, changepoint_signs);
}


BinarySegmentation k_step_bs_wrap(const VecDouble & y, const int n_steps) {
  VecDouble Sy(y.size());
  std::partial_sum(y.begin(), y.end(), Sy.begin());
  return k_step_bs(Sy, n_steps);
}
