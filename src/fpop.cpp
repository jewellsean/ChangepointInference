/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceList.h"
#include <math.h>
#include <stdexcept>
#include <stdlib.h>


void check_min_operation(PiecewiseSquareLoss *cost, PiecewiseSquareLoss *min_prev_cost, PiecewiseSquareLoss *cost_prev, double penalty, int data_i) {
  int status = cost -> check_min_of(min_prev_cost, cost_prev);
  try {
    if(status){
      printf("Lambda = %.20e\n", penalty);
      printf("Error at data_i = %d status = %d\n", data_i, status);
      cost -> set_to_min_env_of(min_prev_cost, cost_prev, true);
      printf("=min_prev_cost\n");
      min_prev_cost -> print();
      printf("=cost_prev + %f\n", penalty);
      cost_prev -> print();
      printf("=new cost model\n");
      cost -> print();
      throw std::runtime_error("min(f, g) error. Please report!");
    }
  } catch(int e) {
    printf("An exception occured %d \n", e);
    throw std::runtime_error("min(f, g) error. Please report!");
  }
}

PiecewiseSquareLoss * fpop_update_i(PiecewiseSquareLoss *cost, PiecewiseSquareLoss *cost_prev, double data, double penalty, int data_i, int verbose) {

  PiecewiseSquareLoss min_prev_cost;

  min_prev_cost.set_to_unconstrained_min_of(cost_prev, verbose);
  min_prev_cost.set_prev_seg_end(data_i - 1);
  min_prev_cost.add(0, 0, penalty);
  cost -> set_to_min_env_of(&min_prev_cost, cost_prev, verbose);
  check_min_operation(cost, &min_prev_cost, cost_prev, penalty, data_i);
  cost -> add(0.5, - data, data * data / 2);
  return(cost);
}


PiecewiseSquareLosses fpop
        (double *data_vec, int data_count, double penalty, double min_mean, double max_mean) {

  PiecewiseSquareLoss *cost, *cost_prev;
  PiecewiseSquareLoss min_prev_cost;
  PiecewiseSquareLosses cost_model_mat(data_count);

  int verbose=0;
  for(int data_i=0; data_i< data_count; data_i++){
    cost = &cost_model_mat[data_i];
    if(data_i==0){
      cost -> piece_list.emplace_back(0.5, - data_vec[0], data_vec[0] * data_vec[0] / 2, min_mean, max_mean, -1, false);
      cost_prev = cost;
    }else {
      cost_prev = fpop_update_i(&cost_model_mat[data_i], cost_prev, data_vec[data_i], penalty, data_i, verbose);
    }
  }
  return cost_model_mat;
}

PiecewiseSquareLosses fpop_custom
        (double *data_vec, int data_count,
         PiecewiseSquareLoss *start_cost,
         double penalty,
         int verbose
        ){

  PiecewiseSquareLoss *cost_prev;
  PiecewiseSquareLosses cost_model_mat(data_count);

  cost_prev = start_cost;

  // build cost functions for each data point, starting with start_cost cost function
  for(int data_i=0; data_i < data_count; data_i++){
    cost_prev = fpop_update_i(&cost_model_mat[data_i], cost_prev, data_vec[data_i], penalty, data_i, verbose);
  }
  return cost_model_mat;
}

void decode_fpop(PiecewiseSquareLosses cost_model_mat,
                 int data_count,
                 double *cost_mat,
                 int *end_vec,
                 double *mean_vec,
                 int *intervals_mat) {

  PiecewiseSquareLoss *cost;
  double best_cost, best_mean, prev_mean;
  int prev_seg_end=data_count;

  for(int i=0; i< data_count; i++){
    cost = &cost_model_mat[i];

/*    printf("cost at data_i = %d\n", i);
    cost -> print();*/

    intervals_mat[i] = cost->piece_list.size();
    cost->Minimize
            (&best_cost, &best_mean,
             &prev_seg_end, &prev_mean);

    cost_mat[i] = best_cost;
  }

  // first step
  cost = &cost_model_mat[data_count - 1];
  cost->Minimize
          (&best_cost, &best_mean,
           &prev_seg_end, &prev_mean);


  int prev_seg_old = data_count - 1;
  int out_i=0;

  // loop over all prev. changepoints
  while(prev_seg_old >= 0){
    if (prev_seg_old < data_count - 1) {
      cost = &cost_model_mat[prev_seg_end];
      cost->Minimize
              (&best_cost, &best_mean,
               &prev_seg_end, &prev_mean);
    }
    for (int t = prev_seg_old; t > prev_seg_end; t--){
      mean_vec[out_i] = best_mean;
      end_vec[out_i] = prev_seg_end;
      out_i++;
    }
    prev_seg_old = prev_seg_end;
  }
}

std::list<int> extract_changepoints(PiecewiseSquareLosses cost_model_mat, int data_count) {
  int prev_seg_old = data_count - 1;
  std::list<int> changepoints;
  // loop over all prev. changepoints
  while(prev_seg_old >= 0){
    prev_seg_old = cost_model_mat[prev_seg_old].getMinIntervalInd();
    changepoints.emplace_front(prev_seg_old);
  }
  return changepoints;
}


PiecewiseBiSquareLosses fpop_2d_custom_start(double *data_vec, int data_count,
                                             PiecewiseSquareLoss *start_cost,
                                             double penalty,
                                             bool forward,
                                             double nuTy,
                                             double norm_constant,
                                             int verbose) {

  PiecewiseBiSquareLosses collection;
  PiecewiseBiSquareLoss start;
  PiecewiseBiSquareLoss min_over_mu;

  start.set_to_pw_u(start_cost);
  collection.piece_list.clear();
  collection.piece_list.emplace_back(start);

  // build cost functions for each data point, starting with start_cost cost function
  for(int data_i=0; data_i < data_count; data_i++){
    min_over_mu = collection.min_u();

    if (verbose) {
      printf("--------\n");
      printf("Collection at data_i=%d\n", data_i);
      printf("Min over mu of all elements in collection at t-1\n");
      min_over_mu.print();
    }

    min_over_mu.add(0, 0, 0, 0, 0, penalty);
    collection.collect(&min_over_mu);

    if (forward) {
      collection.add(0.5,
                     nuTy / norm_constant - data_vec[data_i],
                     - 1 / norm_constant,
                     1 / (2 * norm_constant * norm_constant),
                     (1 / norm_constant) * (data_vec[data_i] - nuTy / norm_constant),
                     (nuTy / norm_constant) * (nuTy / (2 * norm_constant) - data_vec[data_i]) +
                     data_vec[data_i] * data_vec[data_i] / 2);
    } else {
      collection.add(0.5,
                     -nuTy / norm_constant - data_vec[data_i],
                     1 / norm_constant,
                     1 / (2 * norm_constant * norm_constant),
                     (1 / norm_constant) * (-data_vec[data_i] - nuTy / norm_constant),
                     (nuTy / norm_constant) * (nuTy / (2 * norm_constant) + data_vec[data_i]) +
                     data_vec[data_i] * data_vec[data_i] / 2);
    }

    if (verbose) {
      collection.print();
    }

  }
  return collection;
}
