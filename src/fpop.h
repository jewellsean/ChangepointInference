#include "funPieceList.h"
#include <stdexcept>

PiecewiseSquareLosses fpop(double *, int, double, double, double);

PiecewiseSquareLosses fpop_custom(double *, int, PiecewiseSquareLoss *, double, int);

void decode_fpop(PiecewiseSquareLosses, int, double *, int *, double *, int *);

void check_min_operation(PiecewiseSquareLoss *, PiecewiseSquareLoss *, PiecewiseSquareLoss *, double, int);

PiecewiseSquareLoss * fpop_update_i(PiecewiseSquareLoss *, PiecewiseSquareLoss *, double, double, int, int);

PiecewiseBiSquareLosses fpop_2d_custom_start(double *data_vec, int data_count,
                                             PiecewiseSquareLoss *start_cost,
                                             double penalty,
                                             bool forward,
                                             double nuTy,
                                             double norm_constant,
                                             int verbose);

std::list<int> extract_changepoints(PiecewiseSquareLosses cost_model_mat, int data_count);