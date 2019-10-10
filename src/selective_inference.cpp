#include "selective_inference.h"
#include "fpop.h"
#include "funPieceList.h"
#include <vector>
#include "utils.h"
#include <algorithm>    // std::min
#include <math.h>
#include <stdexcept>
/*
 * Contrast vector v
 *
 * v_t = 0,                     if t < t_L or t > t_R
 * v_t = 1 / (thj - t_L + 1),   if t_L <= t <= thj
 * v_t = 1 / (t_R - thj),       if thj + 1 <= t <= t_R
 *
 * for some changepoint of interest thj and endpoints
 *
 * t_L = max(0, thj - window_size + 1)
 * t_R = min(data_count, thj + window_size)
 *
 * and window_size
 *
 */
double construct_vTy(double * sub_data_f, int n_sub_f, double * sub_data_r, int n_sub_r) {
  double vTy = 0;
  double L_norm = 1.0 / n_sub_f;
  double R_norm = 1.0 / n_sub_r;

  for (int i = 0; i < n_sub_f; i++) {
    vTy += sub_data_f[i] * L_norm;
  }

  for (int i = 0; i < n_sub_r; i++) {
    vTy += -sub_data_r[i] * R_norm;
  }

  return vTy;
}

double construct_left_shift(double phi, int n_sub_f, int n_sub_r, double vTy) {
  double D = 1 + (double) n_sub_f / n_sub_r;
  double left_shift = (phi - vTy) / D;
  return left_shift;
}


double construct_right_shift(double phi, int n_sub_f, int n_sub_r, double vTy) {
  double D = 1 + (double) n_sub_r / n_sub_f;
  double right_shift = -(phi - vTy) / D;
  return right_shift;
}


int thj_in_model_at_phi(
        PiecewiseSquareLosses *cost_fwd,
        PiecewiseSquareLosses *cost_rev,
        int thj, // changepoint of interest
        int window_size, // size of window around thj
        int data_count, // number of data points
        double *data_vec, // original data
        double *data_vec_rev, // flipped data
        double phi, // shift to data
        double penalty, // tuning parameter to penalize the number of spikes
        double *opt_cost,
        int verbose
        ) {

  if (window_size < 2) throw;


  int sub_f_start = std::max(thj - window_size + 1, 0);
  int n_sub_f = thj - sub_f_start + 1;
  double * sub_data_f = subset_array(data_vec, sub_f_start, thj + 1);

  int sub_r_start = std::max(data_count - 1 - thj - window_size, 0);
  int n_sub_r = data_count - 1 - thj - 1 - sub_r_start + 1;
  double * sub_data_r = subset_array(data_vec_rev, sub_r_start, data_count - thj - 1);
  double vTy = construct_vTy(sub_data_f, n_sub_f, sub_data_r, n_sub_r);

  double left_shift = construct_left_shift(phi, n_sub_f, n_sub_r, vTy);
  double right_shift = construct_right_shift(phi, n_sub_f, n_sub_r, vTy);
  double * sub_phi_f = add_constant_array(sub_data_f, left_shift, n_sub_f);
  double * sub_phi_r = add_constant_array(sub_data_r, right_shift, n_sub_r);

  // if (verbose) debug_pm_theta(data_vec, data_vec_rev, sub_phi_f, sub_phi_r, data_count, n_sub_f, n_sub_r);

  int cost_f_start = (int) std::max(thj - window_size, 0);
  int cost_r_start = (int) std::max(data_count - 1 - thj - window_size - 1, 0);


  PiecewiseSquareLoss *cost_f_piece, *cost_r_piece;
  PiecewiseSquareLoss fwd_0, rev_0;

  if (thj - window_size < 0) {
    fwd_0.piece_list.emplace_back(0, 0, 0, -INFINITY, INFINITY, 0, 0); // degenerate cost section to update
    cost_f_piece = &fwd_0;

  } else {
    cost_f_piece = &(cost_fwd -> at(cost_f_start));
  }

  if (data_count - 1 - thj - window_size - 1 < 0) {
    rev_0.piece_list.emplace_back(0, 0, 0, -INFINITY, INFINITY, 0, 0);
    cost_r_piece = &rev_0;
  } else {
    cost_r_piece = &(cost_rev -> at(cost_r_start));
  }
  PiecewiseSquareLosses costs_fwd_to_thj = fpop_custom(sub_phi_f, n_sub_f, cost_f_piece, penalty, 0);
  PiecewiseSquareLosses costs_rev_to_thj_1 = fpop_custom(sub_phi_r, n_sub_r, cost_r_piece, penalty, 0);


  // cost of segmenting data with changepoint at thj
  // min. each seperately...if data is random a.s. leads to a changepoint
  double best_cost_no_change_left, best_mean, prev_mean;
  int prev_seg_end;

  if (verbose) {
    printf("Forward cost up to thj \n");
    printf("cost at %d\n", cost_f_start);
    cost_f_piece -> print();

    for (int i = 0; i < n_sub_f; i++) {
      printf("cost at data_i %d\n", i);
      costs_fwd_to_thj[i].print();
    }

  }


  costs_fwd_to_thj[n_sub_f - 1].Minimize(&best_cost_no_change_left, &best_mean, &prev_seg_end, &prev_mean);

  if (verbose) {
    printf("Reverse cost up to thj + 1 \n");
    costs_rev_to_thj_1[n_sub_r - 1].print();
  }


  double best_cost_no_change_right;
  costs_rev_to_thj_1[n_sub_r - 1].Minimize(&best_cost_no_change_right, &best_mean, &prev_seg_end, &prev_mean);

  double best_cost_with_changepoint_at_thj = best_cost_no_change_left + best_cost_no_change_right + penalty;

  // cost of segmenting data with NO changepoint at thj
  // add these two functions together and constrain the mean to be equal
  PiecewiseSquareLoss cost_no_change_thj;
  cost_no_change_thj.set_to_addition_of(&costs_fwd_to_thj[n_sub_f - 1], &costs_rev_to_thj_1[n_sub_r - 1], 0);

  double best_cost_no_change;
  cost_no_change_thj.Minimize(&best_cost_no_change, &best_mean, &prev_seg_end, &prev_mean);

  if (verbose) {
    printf("The best cost with no change at thj \t %f \n", best_cost_no_change);
    printf("The best cost with a change at thj \t %f \n", best_cost_with_changepoint_at_thj);

  }
  if (best_cost_no_change < best_cost_with_changepoint_at_thj) {
    *opt_cost = best_cost_no_change;
    return 0;
  } else {
    *opt_cost = best_cost_with_changepoint_at_thj;
    return 1;
  }
}


void grid_search_single_thj(PiecewiseSquareLosses *cost_fwd,
                 PiecewiseSquareLosses *cost_rev,
                 int thj, // changepoint of interest
                 int window_size, // size of window around thj
                 int data_count, // number of data points
                 double *data_vec, // original data
                 double *data_vec_rev, // flipped data
                 double *phis, // grid of phis to search over
                 int phi_count,  // number of phis to search over
                 int *model_at_phi, // indicator variables denoting whether or not thj is in M(y^phi)
                 double penalty, // tuning parameter to penalize the number of spikes
                 int verbose) {


  double opt_cost;
  for (int phi_i = 0; phi_i < phi_count; phi_i++) {
    int in_out = thj_in_model_at_phi(
            cost_fwd,
            cost_rev,
            thj, // changepoint of interest
            window_size, // size of window around thj
            data_count, // number of data points
            data_vec, // original data
            data_vec_rev, // flipped data
            phis[phi_i], // shift to data
            penalty, // tuning parameter to penalize the number of spikes
            &opt_cost,
            verbose);
    model_at_phi[phi_i] = in_out;
  }
}




PiecewiseSquareLoss thj_in_model(
        PiecewiseSquareLosses *cost_fwd,
        PiecewiseSquareLosses *cost_rev,
        int thj, // changepoint of interest
        int window_size, // size of window around thj
        int data_count, // number of data points
        double *data_vec, // original data
        double *data_vec_rev, // flipped data
        double penalty, // tuning parameter to penalize the number of spikes
        int verbose
) {

  if (window_size < 1) throw;


  int sub_f_start = std::max(thj - window_size + 1, 0);
  int n_sub_f = thj - sub_f_start + 1;
  double * sub_data_f = subset_array(data_vec, sub_f_start, thj + 1);

  int sub_r_start = std::max(data_count - 1 - thj - window_size, 0);
  int n_sub_r = data_count - 1 - thj - 1 - sub_r_start + 1;
  double * sub_data_r = subset_array(data_vec_rev, sub_r_start, data_count - thj - 1);
  double vTy = construct_vTy(sub_data_f, n_sub_f, sub_data_r, n_sub_r);

  double norm_constant_f = 1 + ((double) n_sub_f / n_sub_r);
  double norm_constant_r = 1 + ((double) n_sub_r / n_sub_f);

  int cost_f_start = (int) std::max(thj - window_size, 0);
  int cost_r_start = (int) std::max(data_count - 1 - thj - window_size - 1, 0);

  PiecewiseSquareLoss *cost_f_piece, *cost_r_piece;
  PiecewiseSquareLoss fwd_0, rev_0;

  if (thj - window_size < 0) {
    fwd_0.piece_list.emplace_back(0, 0, 0, MACHINE_MIN, MACHINE_MAX, 0, 0); // degenerate cost section to update
    cost_f_piece = &fwd_0;

  } else {
    cost_f_piece = &(cost_fwd -> at(cost_f_start));
  }

  if (data_count - 1 - thj - window_size - 1 < 0) {
    rev_0.piece_list.emplace_back(0, 0, 0, MACHINE_MIN, MACHINE_MAX, 0, 0);
    cost_r_piece = &rev_0;
  } else {
    cost_r_piece = &(cost_rev -> at(cost_r_start));
  }

  PiecewiseBiSquareLosses fwd_2d = fpop_2d_custom_start(sub_data_f, n_sub_f,cost_f_piece, penalty, 1, vTy, norm_constant_f, verbose);
  PiecewiseBiSquareLosses rev_2d = fpop_2d_custom_start(sub_data_r, n_sub_r,cost_r_piece, penalty, 0, vTy, norm_constant_r, verbose);

  PiecewiseSquareLoss fwd_min, rev_min, c_change_at_thj;
  fwd_min = fwd_2d.min_u().get_univariate_p();
  rev_min = rev_2d.min_u().get_univariate_p();
  c_change_at_thj.set_to_addition_of(&fwd_min, &rev_min, 0);
  c_change_at_thj.add(0, 0, penalty);
  c_change_at_thj.set_prev_seg_end(1); // there is a changepoint at thj

//  printf("eval at baseline (change)= %f\n", c_change_at_thj.findCost(vTy));

  PiecewiseBiSquareLosses c_no_change_at_thj;
  c_no_change_at_thj.set_to_addition_of(&fwd_2d, &rev_2d, 0);

//  printf("addition of fwd and reverse costs\n");
//  c_no_change_at_thj.print();

  PiecewiseSquareLoss c_no_change;
  c_no_change = c_no_change_at_thj.min_u().get_univariate_p();
//  printf("printing c_no_change...min over u\n");
//  c_no_change.print();
//  printf("eval at baseline (no change)= %f\n", c_no_change.findCost(-1.0));

//  PiecewiseBiSquareLoss c_no_change;
//  c_no_change = c_no_change_at_thj.min_u();
//  c_no_change.print();

  c_no_change.set_prev_seg_end(0); // there is no changepoint at thj

  PiecewiseSquareLoss optimal_cost_in_phi;
  optimal_cost_in_phi.set_to_min_env_of(&c_change_at_thj, &c_no_change, 0);

  return optimal_cost_in_phi;

}


void check_selective_inference(PiecewiseSquareLoss * analytic_phi,
        int thj, // changepoint of interest
        int window_size, // size of window around thj
        int data_count, // number of data points
        double *data_vec, // original data
        double penalty, // tuning parameter to penalize the number of spikes
        int verbose) {

  const double MIN = -100;
  const double MAX = 100;

  SquareLossPieceList::iterator it;
  double phi_eval, analytic_cost, manual_cost;

  for (it = analytic_phi->piece_list.begin(); it != analytic_phi->piece_list.end(); it++) {
    phi_eval = MidMean(it -> min_mean, it -> max_mean);

    if (phi_eval > MIN && phi_eval < MAX) {
      analytic_cost = it -> getCost(phi_eval);

      // run fpop on yphi

      int sub_f_start = std::max(thj - window_size + 1, 0);
      int n_sub_f = thj - sub_f_start + 1;
      int sub_r_end = std::min(data_count - 1, thj + window_size);
      int n_sub_r = sub_r_end - thj;

      double norm_constant_f = 1 + ((double) n_sub_f / n_sub_r);
      double norm_constant_r = 1 + ((double) n_sub_r / n_sub_f);

      double vTy = 0;
      for (int i = 0; i < data_count; i++) {
        if (i >= sub_f_start && i <= thj) {
          vTy += data_vec[i] / n_sub_f;
        } if (i > thj && i <= sub_r_end) {
          vTy -= data_vec[i] / n_sub_r;
        }
      }

      PiecewiseSquareLoss *cost_prev;
      PiecewiseSquareLosses cost_model_mat(data_count);
      PiecewiseSquareLoss start;
      start.piece_list.emplace_back(0, 0, 0, -INFINITY, INFINITY, 0, 0); // degenerate cost section to update
      cost_prev = &start;

      // build cost functions for each data point, starting with start_cost cost function
      double next_data_point;
      for(int data_i=0; data_i < data_count; data_i++){

        if (data_i < sub_f_start || data_i > sub_r_end) {
          next_data_point = data_vec[data_i];
        } if (data_i >= sub_f_start && data_i <= thj) {
          next_data_point = data_vec[data_i] + (phi_eval - vTy) / norm_constant_f;
        } if (data_i > thj && data_i <= sub_r_end) {
          next_data_point = data_vec[data_i] - (phi_eval - vTy) / norm_constant_r;
        }

        cost_prev = fpop_update_i(&cost_model_mat[data_i], cost_prev, next_data_point, penalty, data_i, verbose);
      }

      manual_cost = cost_model_mat[data_count-1].getMinCost();

      if (ABS(manual_cost - analytic_cost) > DIFF_EPS) {
        printf("analytic cost incorrect. different between analytic cost and manual cost at phi (%f)= \t %.50f", phi_eval, analytic_cost - manual_cost);
        throw std::runtime_error("analytic cost incorrect. Please report!");
      }

    }


  }

}


double construct_vTy(int thj, int window_size, double * data_vec, int data_count) {
  int sub_f_start = std::max(thj - window_size + 1, 0);
  int n_sub_f = thj - sub_f_start + 1;
  int sub_r_end = std::min(data_count - 1, thj + window_size);
  int n_sub_r = sub_r_end - thj;

  double vTy = 0;
  for (int i = 0; i < data_count; i++) {
    if (i >= sub_f_start && i <= thj) {
      vTy += data_vec[i] / n_sub_f;
    } if (i > thj && i <= sub_r_end) {
      vTy -= data_vec[i] / n_sub_r;
    }
  }

  return vTy;

}

double construct_nu2(int thj, int window_size, int data_count) {
  int sub_f_start = std::max(thj - window_size + 1, 0);
  int n_sub_f = thj - sub_f_start + 1;
  int sub_r_end = std::min(data_count - 1, thj + window_size);
  int n_sub_r = sub_r_end - thj;
  return ((1.0 / n_sub_f) + (1.0 / n_sub_r));
}

double calc_p_value(PiecewiseSquareLoss * analytic_phi,
                    int thj, // changepoint of interest
                    int window_size, // size of window around thj
                    int data_count, // number of data points
                    double *data_vec, // original data
                    double sig, // noise variance
                    int verbose) {



  double vTy = construct_vTy(thj, window_size, data_vec, data_count);
  double nu_norm = construct_nu2(thj, window_size, data_count);
  SquareLossPieceList::iterator it;

  // numerically safe
  double n1 = -INFINITY;
  double d1 = -INFINITY;
  double arg2;
  for (it = analytic_phi->piece_list.begin(); it != analytic_phi->piece_list.end(); it++) {
    if (it->data_i == 1) { // this segment is contained
      double a, b;
      a = pnorm_log(it -> max_mean / sqrt(nu_norm * sig));
      b = pnorm_log(it -> min_mean / sqrt(nu_norm * sig));
      arg2 = log_subtract(a, b);
      d1 = log_sum_exp(d1, arg2);

      if (it->max_mean >= ABS(vTy)) {
        arg2 = log_subtract(pnorm_log(it -> max_mean / sqrt(nu_norm * sig)),
                            pnorm_log(std::max(it -> min_mean, ABS(vTy)) / sqrt(nu_norm * sig)));
        n1 = log_sum_exp(n1, arg2);
      }
      if (it->min_mean <= -1 * ABS(vTy)) {
        arg2 = log_subtract(pnorm_log(std::min(it -> max_mean, -ABS(vTy)) / sqrt(nu_norm * sig)),
                            pnorm_log(it -> min_mean / sqrt(nu_norm * sig)));
        n1 = log_sum_exp(n1, arg2);
      }
    }
  }
  return (exp(n1 - d1));
}