#include "fpop.h"
#include "funPieceList.h"
#include "fpop_inference.h"
#include "selective_inference.h"
#include "utils.h"
#include <vector>
#include <Rcpp.h>
#include <algorithm>    // std::min
#include <math.h>
#include "binary_segmentation.h"
#include "bs_inference.h"

using namespace Rcpp;

/*
 *
 * Utility functions
 *
 */

// utility function for fpop intervals
NumericMatrix convert_loss(PiecewiseSquareLoss *p, int s) {
  const int ncols = 8;
  const int n_pieces = p -> piece_list.size();
  NumericMatrix out(n_pieces, ncols);
  
  int row_i = 0;
  for (SquareLossPieceList::iterator it = p -> piece_list.begin(); it != p-> piece_list.end(); it++) {
    out(row_i, 0) = it -> Square;
    out(row_i, 1) = it -> Linear;
    out(row_i, 2) = it -> Constant;
    out(row_i, 3) = it -> min_mean;
    out(row_i, 4) = it -> max_mean;
    out(row_i, 5) = it -> prev_mean;
    out(row_i, 6) = it -> data_i + 1;
    out(row_i, 7) = s;
    row_i++;
  }
  return out;
}

// utility function for BS intervals
NumericMatrix phi_to_mat(const std::vector<PhiInterval> & phi) {
  NumericMatrix out_mat(phi.size(), 3);
  for (int row_i = 0; row_i < phi.size(); row_i++) {
    out_mat(row_i, 0) = phi[row_i].L;
    out_mat(row_i, 1) = phi[row_i].U;
    out_mat(row_i, 2) = phi[row_i].contained;
  }
  return out_mat;
}

/*
 *
 * Estimation functions
 *
 */


// [[Rcpp::export(name = ".fpop")]]
List fpop_interface2
        (NumericVector data,
         double penalty,
         double min_mean,
         double max_mean,
         NumericVector cost_mat_r,
         IntegerVector end_vec_r,
         NumericVector mean_vec_r,
         IntegerVector intervals_mat_r){


  double *data_ptr = data.begin();
  int data_count = data.size();
  double *cost_mat = cost_mat_r.begin();
  int *end_vec = end_vec_r.begin();
  double *mean_vec = mean_vec_r.begin();
  int *intervals_mat = intervals_mat_r.begin();

  PiecewiseSquareLosses out = fpop(data_ptr, data_count, penalty, min_mean, max_mean);
  decode_fpop(out, data_count, cost_mat, end_vec, mean_vec, intervals_mat);
  
  List pw_losses(out.size());
  for (int i = 0; i < out.size(); i++) {
    pw_losses[i] = convert_loss(&out[i], i + 1);
  }
  
  return pw_losses;

}


// [[Rcpp::export(name = ".k_step_bs")]]
NumericMatrix k_step_bs_r(std::vector<double> y, int n_steps) {
  BinarySegmentation out = k_step_bs_wrap(y, n_steps);
  NumericMatrix out_mat(n_steps, 2);
  for (int row_i = 0; row_i < n_steps; row_i++) {
    out_mat(row_i, 0) = out.ordered_changepoints[row_i] + 1; // indexing difference btwn c++ and R
    out_mat(row_i, 1) = out.changepoint_signs[row_i];
  }
  return out_mat;
}

/*
 *
 * Inference functions
 *
 */

// [[Rcpp::export(name = ".fpop_inference")]]
List fpop_inference_interface_recycle
        (NumericVector data,
         double penalty,
         int window_size,
         double sig,
         int return_dev = 0) {


  double *data_ptr = data.begin();
  int data_count = data.size();
  int verbose = 0;
  double pval;


// forward pass
  PiecewiseSquareLosses cost_model_fwd = fpop(data_ptr, data_count, penalty, MACHINE_MIN, MACHINE_MAX);

// backward pass
  double *data_vec_rev = reverse_data(data_ptr, data_count);
  PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, penalty, MACHINE_MIN, MACHINE_MAX);

  std::list<int> ll = extract_changepoints(cost_model_fwd, data_count);
  std::list<int>::iterator it;


//  for each changepoint determine pval
  const int ncols = 3; // changepoint + 1 (in R notation), pval, approximation error
  const int nrows = ll.size() - 1;
  NumericMatrix out_mat(nrows, ncols);
  List phi_intervals(nrows);

  int row_i = 0;
  for (it = ll.begin(); it != ll.end(); ++it) {

    try {
      if (*it > 0) {
        FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, data_ptr, data_count, penalty, *it, window_size, sig);

        out_mat(row_i, 0) = out.thj + 1;
        out_mat(row_i, 1) = out.pval;
        out_mat(row_i, 2) = out.approximation_error;

        if (return_dev) {
          phi_intervals[row_i] = convert_loss(&out.model, 0);
        } else {
          phi_intervals[row_i] = nullptr;
        }

        row_i++;
      }
    } catch (std::exception &ex) {
      forward_exception_to_r(ex);
    } catch (...) {
      ::Rf_error("c++ exception (unknown reason)");
    }
  }

  return List::create(out_mat, phi_intervals);

}

 // [[Rcpp::export(name = ".k_step_bs_inference")]]
 List k_step_bs_inference_r(std::vector<double> y, int n_steps, int type, double sigma, int window_size = 0, double approximation_threshold = 5, int return_dev = 0) {
   BSInference out = k_step_bs_inference(y, n_steps, type, sigma, window_size, approximation_threshold);
   NumericMatrix out_mat(n_steps, 3);
   List phi_intervals(n_steps);

   for (int row_i = 0; row_i < n_steps; row_i++) {
     out_mat(row_i, 0) = out.sorted_changepoints[row_i + 1] + 1; // indexing difference btwn c++ and R
     out_mat(row_i, 1) = out.pvals[row_i];
     out_mat(row_i, 2) = out.approximation_errors[row_i];

     if (return_dev) {
       phi_intervals[row_i] = phi_to_mat(out.phi_intervals[row_i]);
     } else {
       phi_intervals[row_i] = nullptr;
     }

   }
   return List::create(out_mat, phi_intervals, out.ordered_changepoints);
 }