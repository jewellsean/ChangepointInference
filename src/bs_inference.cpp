#include "binary_segmentation.h"
#include "bs_inference.h"
#include <algorithm>
#include <math.h>
#include <numeric>
#include "utils.h"
// class defns

Contrast::Contrast(int l, int t, int r, const VecDouble & Sy) {
  thj = t;
  tL = l;
  tR = r;
  // contrast aspects to store
  // 1 / vTv and
  // vTy
  double v1 = 1.0 / (thj - tL + 1);
  double v2 = -1.0 / (tR - thj);
  nuTnu_inv = (1.0 / ((v1 * v1) * (thj - tL + 1)  + (v2 * v2) * (tR - thj)));

  if (tL == 0) {
    nuTy = (Sy[thj] - 0) / (thj - tL + 1) - (Sy[tR] - Sy[thj]) / (tR - thj);
  } else
  {
    nuTy = (Sy[thj] - Sy[tL - 1]) / (thj - tL + 1) - (Sy[tR] - Sy[thj]) / (tR - thj);
  }
  shift = nuTy;
}

void Contrast::set_shift(double a) {shift = a;}

Contrast::Contrast() {}

void Contrast::update_nuTY(const VecDouble & y) {
  VecDouble Sy(y.size());
  std::partial_sum(y.begin(), y.end(), Sy.begin());
  if (tL == 0) {
    nuTy = (Sy[thj] - 0.0) / (thj - tL + 1) - (Sy[tR] - Sy[thj]) / (tR - thj);
  } else
  {
    nuTy = (Sy[thj] - Sy[tL - 1]) / (thj - tL + 1) - (Sy[tR] - Sy[thj]) / (tR - thj);
  }
}

IntervalInspection::IntervalInspection(int s, int e, int t, int d) {
  start = s;
  end = e;
  tau = t;
  change_sign = d;
}

PhiInterval::PhiInterval(double l, double u, int c) {
  L = l;
  U = u;
  contained = c;
}

PhiInterval::PhiInterval() {}

void PhiInterval::print() {
  printf("%f \t %f \t %d\n", L, U, contained);
}

BSInferenceOpts BSInferenceOpts_from_simple(int t, int n) {
  return BSInferenceOpts(t, n, -1e3, 1e3, 0.1, 20, 1e-8);
}

BSInferenceOpts::BSInferenceOpts() {}

BSInferenceOpts::BSInferenceOpts(int t, int n, double l, double u, double d, int m, double tt) {
  type = t;
  n_steps = n;
  lower_threshold = l;
  upper_threshold = u;
  delta = d;
  maxit = m;
  TOL = tt;
}

void BSInference::add_changepoint_inference(double pval, double approximation_error, std::vector<PhiInterval> phi_interval) {
  pvals.emplace_back(pval);
  phi_intervals.emplace_back(phi_interval);
  approximation_errors.emplace_back(approximation_error);
}

void BSInference::print(){
  printf("Binary segmentation with %d changepoints\n", static_cast<int>(ordered_changepoints.size()));
  printf("...using type = %d\n", opts.type);
  printf("Sorted changepoints (pval):\n");
  int n_pvals = static_cast<int>(sorted_changepoints.size() - 1);
  for(int i = 1; i < n_pvals; i++) {
    printf("changepoint (pval) = %d \t (%f) \n", sorted_changepoints[i], pvals[i - 1]);
  }
}

// class utility
bool partition_sort(PhiInterval a, PhiInterval b) {
  return (a.L < b.L);
}


// end of class defns

Contrast Contrast_from_window(int t, int window_size, const VecDouble & Sy) {
  int n = Sy.size();
  int tL = std::max(t - window_size + 1, 0);
  int tR = std::min(t + window_size, n - 1);
  return Contrast(tL, t, tR, Sy);
}


// compute matrix multiplication
// g^T_(s, tau, e) %*% nu_(thj)
double mult_nu_g(const IntervalInspection & interval, const Contrast & nu) {
  int tL = nu.tL;
  int tR = nu.tR;
  int thj = nu.thj;
  int s = interval.start;
  int e = interval.end;
  int tau = interval.tau;

  double v1 = 1.0 / (thj - tL + 1);
  double v2 = -1.0 / (tR - thj);

  double norm_constant = sqrt(static_cast<double>((tau + 1 - s) * (e - tau)) / (e - s + 1));
  double g1 = -norm_constant / (tau - s + 1);
  double g2 = norm_constant / (e - tau);

  double m1 = std::max(0, std::min(tau, thj) - std::max(s, tL) + 1);
  double m2 = std::max(0, std::min(tau, tR) - std::max(s, thj + 1) + 1);
  double m3 = std::max(0, std::min(e, thj) - std::max(tau + 1, tL) + 1);
  double m4 = std::max(0, std::min(e, tR) - std::max(tau + 1, thj + 1) + 1);

  return m1 * g1 * v1 + m2 * g1 * v2 + m3 * g2 * v1 + m4 * g2 * v2;
}


VecDouble determine_ind_interval(const VecDouble & Sy,
                            const IntervalInspection & int1,
                            const IntervalInspection & int2,
                            const Contrast & nu,
                            int type) {
  int sign_mult = 1;
  if (type == 0) {
    sign_mult = -1;
  }

  double G1v = mult_nu_g(int1, nu);
  double G2v = mult_nu_g(int2, nu);
  double G1y = calculate_csum(Sy, int1.tau, int1.start, int1.end);
  double G2y = calculate_csum(Sy, int2.tau, int2.start, int2.end);

  double Z = (int1.change_sign * G1v + sign_mult * G2v) * nu.nuTnu_inv;
  int sgn = (Z > 0) - (Z < 0);

  double interval_boundary = (Z * nu.nuTy - (int1.change_sign * G1y + sign_mult * G2y)) / Z;

  VecDouble out_interval(2);
  if (sgn > 0) {
    out_interval[0] = interval_boundary;
    out_interval[1] = INFINITY;
  } else if (sgn < 0){
    out_interval[0] = -INFINITY;
    out_interval[1] = interval_boundary;
  } else {
    out_interval[0] = -INFINITY;
    out_interval[1] = INFINITY;
  }
  if ((out_interval[0] > nu.nuTy) || (out_interval[1] < nu.nuTy)) {
    printf("Shift = %f is not contained in obtained interval [%f, %f]\n", nu.nuTy, out_interval[0], out_interval[1]);
    throw std::runtime_error("Shift not contained in obtained interval");
  }
  return out_interval;
}

VecDouble interval_intersection(VecDouble a, VecDouble b) {
  double a_out = std::max(a[0], b[0]);
  double b_out = std::min(a[1], b[1]);
  if (b_out < a_out) {
    throw std::runtime_error("intervals are not overlapping!");
  }
  return VecDouble({a_out, b_out});
}

VecDouble determine_interval(const VecDouble & Sy,
                                  const IntervalInspection & int1,
                                  const IntervalInspection & int2,
                                  const Contrast & nu) {
  VecDouble pos_int = determine_ind_interval(Sy, int1, int2, nu, 0);
  VecDouble neg_int = determine_ind_interval(Sy, int1, int2, nu, 1);
  return interval_intersection(pos_int, neg_int);
}

VecDouble interval_over_region(const VecDouble & Sy,
                               const IntervalInspection & int1,
                               IntervalInspection int2,
                               const Contrast & nu) {

  VecDouble base_interval = VecDouble({-INFINITY, INFINITY});
  VecDouble int_tmp;

  for (int tj = int2.start; tj < int2.end; tj++) {
    if (tj != int1.tau) {
      int2.tau = tj;
      int_tmp = determine_interval(Sy, int1, int2, nu);
//      if (int_tmp[0] > -INFINITY || int_tmp[1] < INFINITY) {
//        printf("\t \t ** tj = %d \t [%f, %f]\n", tj, int_tmp[0], int_tmp[1]);
//      }
      base_interval = interval_intersection(int_tmp, base_interval);
    }
  }

  return base_interval;
}

// modified from https://stackoverflow.com/questions/469477/find-nearest-points-in-a-vector
std::pair<int, int> get_nearest_changepoints(const VecInt & vec, const int & val) {

  auto it = std::lower_bound(vec.begin(), vec.end(), val);

  if (it == vec.end())
    return std::make_pair(vec.back(), vec.back());
  else if (it == vec.begin())
    return std::make_pair(vec.front(), vec.front());
  else
    return std::make_pair(*(it - 1), *(it));
}



VecDouble interval_over_segmentation(const BinarySegmentation & segmentation, const VecDouble & Sy, const Contrast & nu) {
  int n_steps = segmentation.ordered_changepoints.size();
  VecDouble int_region;
  for (int step_i = 0; step_i < n_steps; step_i++) {
//    printf("------\n");
//    printf("at step i = %d\n", step_i);



    if (step_i == 0) {
      IntervalInspection int1 = IntervalInspection(0, Sy.size() - 1,
                                                   segmentation.ordered_changepoints[0], segmentation.changepoint_signs[0]);
      IntervalInspection int2 = IntervalInspection(0, Sy.size() - 1,
                                                   segmentation.ordered_changepoints[0], segmentation.changepoint_signs[0]);
      int_region = interval_over_region(Sy, int1, int2, nu);
    } else {
      int thj = segmentation.ordered_changepoints[step_i];
      int offset = step_i;
      auto first = segmentation.ordered_changepoints.begin();
      auto last = segmentation.ordered_changepoints.begin() + offset;
      VecInt sorted_changepoints(first, last);
      sorted_changepoints.emplace_back(-1);
      sorted_changepoints.emplace_back(Sy.size() - 1);
      std::sort(sorted_changepoints.begin(), sorted_changepoints.end());

      std::pair<int, int> s_e_thj = get_nearest_changepoints(sorted_changepoints, thj);

      int s_thj = s_e_thj.first + 1;
      int e_thj = s_e_thj.second;

      if (e_thj > s_thj) {
        double csum_at_thj = calculate_csum(Sy, thj, s_thj, e_thj);
        int d = (csum_at_thj > 0) - (csum_at_thj < 0);

        IntervalInspection int1 = IntervalInspection(s_thj, e_thj, thj, d);


        for(int seg_i = 0; seg_i <= step_i; seg_i++) {
          int s = sorted_changepoints[seg_i] + 1;
          int e = sorted_changepoints[seg_i + 1];
          if (s < e) {
            IntervalInspection int2 = IntervalInspection(s, e, thj, 0);
            VecDouble int_tmp = interval_over_region(Sy, int1, int2, nu);
            int_region = interval_intersection(int_region, int_tmp);
//            printf("\t seg_i = %d \t [%f, %f]\n", seg_i, int_tmp[0], int_tmp[1]);
          }
        }
      }
    }
//    printf("[%f, %f]\n", int_region[0], int_region[1]);
  }
  return int_region;
}

int check_contained(const BinarySegmentation & segmentation, const BinarySegmentation & orig_segmentation,
        const BSInferenceOpts & opts, const Contrast & nu) {

  if (opts.type == 0) { // thj \in M(y^phi)
    for (auto x : segmentation.ordered_changepoints) {
      if (x == nu.thj) {
        return 1;
      }
    }
    return 0;
  } else if (opts.type == 1) { // M(y) \in M(y^phi) where M encodes the order and changepoint signs

    if (segmentation.ordered_changepoints == orig_segmentation.ordered_changepoints &&
    segmentation.changepoint_signs == orig_segmentation.changepoint_signs) {
      return 1;
    } else {
      return 0;
    }

  } else if (opts.type == 2) { // M(y) \in M(y^phi) where M encodes the order

    if (segmentation.ordered_changepoints == orig_segmentation.ordered_changepoints) {
      return 1;
    } else {
      return 0;
    }

  } else if (opts.type == 3) { // M(y) \in M(y^phi) where M encodes the changepoints

//    printf("====== checking contained function ======= \n");
//    printf("original segmentation \n");
//    for (auto x : orig_segmentation.ordered_changepoints) {printf("%d \t", x);}
//    printf("\n new segmentation \n");
//    for (auto x : segmentation.ordered_changepoints) {printf("%d \t", x);}
//    printf("\n");

    for (auto x : orig_segmentation.ordered_changepoints) {
      if(std::find(segmentation.ordered_changepoints.begin(), segmentation.ordered_changepoints.end(), x) !=
                segmentation.ordered_changepoints.end()) {
      } else {
//        printf("NOT contained\n ----\n");
        return 0;
      }
    }
//    printf("contained\n----\n");
    return 1;
  } else {
    throw std::range_error("type not implemented");
  }

}


PhiInterval
segment_model(const VecDouble &y, const Contrast &nu, const BSInferenceOpts & opts, const BinarySegmentation & orig_segmentation) {
  VecDouble Sy(y.size());
  std::partial_sum(y.begin(), y.end(), Sy.begin());
  BinarySegmentation segmentation = k_step_bs(Sy, opts.n_steps);
  VecDouble interval = interval_over_segmentation(segmentation, Sy, nu);

  //  printf("**** new segment model at phi = %f ***\n", nu.shift);
  //  segmentation.print();

  int contained = check_contained(segmentation, orig_segmentation, opts, nu);
  return PhiInterval(interval[0], interval[1], contained);
}

VecDouble construct_yphi(const VecDouble & y, const Contrast & nu, double shift) {
  VecDouble yphi(y); // deep copy
  double d1 = 1.0 + double (nu.thj - nu.tL + 1) / (nu.tR - nu.thj);
  double d2 = 1.0 + double (nu.tR - nu.thj) / (nu.thj - nu.tL + 1);
  for (int i = nu.tL; i <= nu.thj; i++) {
    yphi[i] = yphi[i] + shift / d1;
  }
  for (int i = nu.thj + 1; i <= nu.tR; i++) {
    yphi[i] = yphi[i] - shift / d2;
  }
  return yphi;
}

void check_phi_interval(double phi_shift, PhiInterval cur_phi_interval) {
  if((phi_shift > cur_phi_interval.U) || (phi_shift < cur_phi_interval.L)) {
    printf("****=======***\n");
    printf("phi_shift = %f\n", phi_shift);
    printf("versus cur interval\n");
    cur_phi_interval.print();
    throw std::runtime_error("phi not bounded\n");
  }
}


std::vector<PhiInterval> partition_phi_space(const VecDouble & y, Contrast & nu, const BSInferenceOpts & opts,
        const BinarySegmentation & orig_segmentation) {
  int verbose = 0;
  double TOL = opts.TOL;
  double lower_thresh = opts.lower_threshold;
  double upper_thresh = opts.upper_threshold;
  double delta = opts.delta;
  int MAXIT = opts.maxit;

  const double nuTy_orig_data = nu.nuTy;
  std::vector<PhiInterval> intervals;

  intervals.emplace_back(segment_model(y, nu, opts, orig_segmentation));

  if (opts.type == 1) {
    return intervals;
  }

  double Lcur = intervals[0].L;
  double Ucur = intervals[0].U;
  VecDouble yphi;
  PhiInterval cur_phi_interval;
  int it = 0;

  while((Lcur > lower_thresh)  || (Ucur < upper_thresh)) {
    if (Lcur > lower_thresh) {
      yphi = construct_yphi(y, nu, Lcur - delta - nuTy_orig_data);
      nu.set_shift(Lcur - delta - nuTy_orig_data);
      nu.update_nuTY(yphi);
      cur_phi_interval = segment_model(yphi, nu, opts, orig_segmentation);

      check_phi_interval(nu.nuTy, cur_phi_interval);

      if (verbose) {
        printf("Segmenting model....\n Current step \n");
        cur_phi_interval.print();
      }

      it = 0;
    while((ABS(cur_phi_interval.U - Lcur) > TOL) && (it < MAXIT)) {
      delta = delta / 2;
      double shift = Lcur - delta;
      nu.set_shift(shift - nuTy_orig_data);
      yphi = construct_yphi(y, nu, shift - nuTy_orig_data);
      nu.update_nuTY(yphi);
      cur_phi_interval = segment_model(yphi, nu, opts, orig_segmentation);
      check_phi_interval(nu.nuTy, cur_phi_interval);
      it = it + 1;
    }
    delta = opts.delta;
    intervals.emplace_back(cur_phi_interval);
      if (verbose) {
        printf("Segmenting model....\n Current step \n");
        cur_phi_interval.print();
      }
    Lcur = cur_phi_interval.L;
  }
  if (Ucur < upper_thresh) {
    yphi = construct_yphi(y, nu, Ucur + delta - nuTy_orig_data);
    nu.set_shift(Ucur + delta - nuTy_orig_data);
    nu.update_nuTY(yphi);
    cur_phi_interval = segment_model(yphi, nu, opts, orig_segmentation);
    check_phi_interval(nu.nuTy, cur_phi_interval);

    it = 0;
    while((ABS(cur_phi_interval.L - Ucur) > TOL) && (it < MAXIT)) {
      delta = delta / 2;
      double shift = Ucur + delta;
      yphi = construct_yphi(y, nu, shift - nuTy_orig_data);
      nu.set_shift(shift - nuTy_orig_data);
      nu.update_nuTY(yphi);
      cur_phi_interval = segment_model(yphi, nu, opts, orig_segmentation);
      check_phi_interval(nu.nuTy, cur_phi_interval);
      it = it + 1;
    }
    delta = opts.delta;
    intervals.emplace_back(cur_phi_interval);
    Ucur = cur_phi_interval.U;
  }
  }
  return intervals;
}

std::vector<PhiInterval> consolidate_phi_intervals(const std::vector<PhiInterval> & intervals, double TOL) {
  if (intervals.size() == 1) {
    return intervals;
  }

  double curL = intervals[0].L;
  int curContained = intervals[0].contained;
  std::vector<PhiInterval> out;

  for (int i = 1; i < static_cast<int>(intervals.size()); i++) {

    if (ABS(intervals[i].L - intervals[i - 1].U) > TOL) {
      return intervals; // not contig.!
    }

    if ((intervals[i].contained != curContained)  && (i < static_cast<int>(intervals.size() - 1))) {
      out.emplace_back(PhiInterval(curL, intervals[i].L, curContained));
      curL = intervals[i].L;
      curContained = intervals[i].contained;
    } else if (i == static_cast<int>(intervals.size() - 1)) {
      if (intervals[i].contained != curContained) {
        out.emplace_back(PhiInterval(curL, intervals[i].L, curContained));
        out.emplace_back(intervals[i]);
      } else {
        out.emplace_back(PhiInterval(curL, intervals[i].U, curContained));
      }
    }
  }
  return out;
}

double calc_prob_set_log(const std::vector<PhiInterval> & intervals,
                     const Contrast nu,
                     double sig, // noise variance
                     int verbose) {
  double vTy = nu.nuTy;
  double nu_norm = 1.0 / nu.nuTnu_inv;

  // numerically safe
  double d1 = -INFINITY;
  double arg2;
  for (auto x : intervals) {
    if (x.contained == 1) {
      double a, b;
      a = pnorm_log(x.U / sqrt(nu_norm * sig));
      b = pnorm_log(x.L / sqrt(nu_norm * sig));
      arg2 = log_subtract(a, b);
      d1 = log_sum_exp(d1, arg2);
    }
  }
  return d1;
}


double calc_p_value(const std::vector<PhiInterval> & intervals,
                    const Contrast nu,
                    double sig, // noise variance
                    int verbose) {
  double vTy = nu.nuTy;
  double nu_norm = 1.0 / nu.nuTnu_inv;

//  printf("vTy = %f \t vTy / var = %f \n", vTy, vTy / sqrt(nu_norm * sig));
  // if number of intervals == 1
  // tail bound to prevent numerical instability
  const double tail_threshold = 30.0;
  if (intervals.size() == 1 && ABS(vTy) >= tail_threshold * sqrt(nu_norm * sig)) {
      double a = 0;
      PhiInterval region = intervals.at(0);
      if (vTy < 0 && region.U < 0) {
        a = ABS(region.U) / sqrt(nu_norm * sig);
      }
      if (vTy > 0 && region.L > 0){
        a = region.L / sqrt(nu_norm * sig);
      }

      if (a > 0) {
        double vDnorm = std::max(ABS(vTy) / sqrt(nu_norm * sig), a);
        return std::min((exp((a * a - (vDnorm * vDnorm)) / 2.0 ) * (a * a + 1)) / (a * vDnorm), 1.0);
      }

  }


  double n1 = -INFINITY;
  double d1 = -INFINITY;
  double arg2;
  for (auto x : intervals) {
    if (x.contained == 1) {
      double a, b;
      a = pnorm_log(x.U / sqrt(nu_norm * sig));
      b = pnorm_log(x.L / sqrt(nu_norm * sig));
      arg2 = log_subtract(a, b);
      d1 = log_sum_exp(d1, arg2);
      if (x.U >= ABS(vTy)) {
        arg2 = log_subtract(pnorm_log(x.U / sqrt(nu_norm * sig)),
                            pnorm_log(std::max(x.L, ABS(vTy)) / sqrt(nu_norm * sig)));
        n1 = log_sum_exp(n1, arg2);
      }
      if (x.L <= -1 * ABS(vTy)) {
        arg2 = log_subtract(pnorm_log(std::min(x.U, -ABS(vTy)) / sqrt(nu_norm * sig)),
                            pnorm_log(x.L / sqrt(nu_norm * sig)));
        n1 = log_sum_exp(n1, arg2);
      }
    }
  }
  return (exp(n1 - d1));
}

void check_regions(const std::vector<PhiInterval> & intervals, const Contrast & nu, const VecDouble & y,
        BinarySegmentation & orig_segmentation, const BSInferenceOpts & opts) {
  double phi_to_check;
  for (auto x : intervals) {

    if (x.L > -INFINITY && x.U < INFINITY) {

      double lower = std::max(opts.lower_threshold, x.L);
      double upper = std::min(opts.upper_threshold, x.U);
      if (x.L > -INFINITY && x.U < INFINITY) {
        phi_to_check = (upper + lower) / 2;
      } else if (x.L > -INFINITY) {
        phi_to_check = lower + 1; // x.U == INF
      } else if (x.U < INFINITY) {
        phi_to_check = upper - 1; // x.L == -INF
      } else {
        phi_to_check = 0;
      }

      if (phi_to_check < x.L || phi_to_check > x.U) { throw std::runtime_error("bad range"); }

      VecDouble yphi = construct_yphi(y, nu, phi_to_check - nu.nuTy);
      BinarySegmentation out = k_step_bs_wrap(yphi, opts.n_steps);

      if (x.contained != check_contained(out, orig_segmentation, opts, nu)) {
        printf("========== \n");
        printf("Error in segmented model for model type %d\n", opts.type);
        if (opts.type == 0) {
          printf("at changepoint = %d\n", nu.thj);
        }
        printf("Checking at phi = %f\n", phi_to_check);
        printf("L = %f \t U = %f \t contained = %d\n", x.L, x.U, x.contained);
        printf("versus contained check =  %d\n", check_contained(out, orig_segmentation, opts, nu));
        printf("original model = \n");
        orig_segmentation.print();
        printf("model at y_phi segmentation\n");
        out.print();
        throw std::runtime_error("failed check\n");
      }

    }
  }
}

BSInference k_step_bs_inference(const VecDouble & y, const int n_steps, const int type, const double sigma, const int window_size = 0, const double threshold = 5) {
  int verbose = 0;
  VecDouble Sy(y.size());
  std::partial_sum(y.begin(), y.end(), Sy.begin());
  BinarySegmentation segmentation = k_step_bs(Sy, n_steps);
  BSInferenceOpts opts = BSInferenceOpts_from_simple(type, n_steps);
  BSInference out = BSInference(segmentation.ordered_changepoints, segmentation.changepoint_signs, opts, y.size());

  Contrast nu;
  for (int i = 1; i <= n_steps; i++) {
    if (type == 0) {
       nu = Contrast_from_window(out.sorted_changepoints[i], window_size, Sy);
    } else {
      nu = Contrast(out.sorted_changepoints[i - 1] + 1, out.sorted_changepoints[i], out.sorted_changepoints[i + 1], Sy);
    }

    if (verbose) {
      printf("====================\n\n\n\n");
      printf("checking changepoint = %d\n", out.sorted_changepoints[i]);
    }
    double vTy = nu.nuTy;
    double nu_norm = 1.0 / nu.nuTnu_inv;
    double lower = std::min(- threshold * sqrt(nu_norm * sigma), -ABS(vTy));
    double upper = std::max(threshold * sqrt(nu_norm * sigma), ABS(vTy));

    BSInferenceOpts opts = BSInferenceOpts(type, n_steps, lower, upper, 0.1, 20, 1e-8);
    std::vector<PhiInterval> partition_out = partition_phi_space(y, nu, opts, out);
    std::sort(partition_out.begin(), partition_out.end(), partition_sort);
    std::vector<PhiInterval> out_consolidated = consolidate_phi_intervals(partition_out, opts.TOL);
    nu.update_nuTY(y);

    double pr_bprime = calc_prob_set_log(out_consolidated, nu, sigma, 0);
    double pr_xtra_regions = -INFINITY;

    if (opts.type != 1) {
      // add conservative p-values here
      //    basically need a partition from -inf, min_value, 1
      if (out_consolidated.at(0).L > -INFINITY) {
        PhiInterval lower_phi = PhiInterval(-INFINITY, out_consolidated.at(0).L, 1);
        out_consolidated.insert(out_consolidated.begin(), lower_phi);

        std::vector<PhiInterval> extra_partition;
        extra_partition.emplace_back(lower_phi);
        pr_xtra_regions = log_sum_exp(pr_xtra_regions, calc_prob_set_log(extra_partition, nu, sigma, 0));
      }

      if (out_consolidated.at(out_consolidated.size() - 1).U < INFINITY) {
        //    basically need a partition from max_value, inf, 1
        PhiInterval upper_phi = PhiInterval(out_consolidated.at(out_consolidated.size() - 1).U, INFINITY, 1);
        out_consolidated.emplace_back(upper_phi);

        std::vector<PhiInterval> extra_partition;
        extra_partition.emplace_back(upper_phi);
        pr_xtra_regions = log_sum_exp(pr_xtra_regions, calc_prob_set_log(extra_partition, nu, sigma, 0));
      }
    }

    double approximation_error = exp(log_subtract(pr_xtra_regions, log_sum_exp(pr_bprime, pr_xtra_regions)));

    check_regions(out_consolidated, nu, y, out, opts);
    double pval = calc_p_value(out_consolidated, nu, sigma, 0);

    if (verbose) {
      printf("pval (error) = %f \t (%f) \t thj = %d \n", pval, approximation_error, out.sorted_changepoints[i]);
      printf("consolidated regions\n");
      for (auto x : out_consolidated) {
        x.print();
      }
    }

    out.add_changepoint_inference(pval, approximation_error, out_consolidated);
  }
  return out;
}