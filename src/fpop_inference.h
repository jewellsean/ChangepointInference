#include "funPieceList.h"
#include <stdexcept>

class FpopInference {
public:
    double pval;
    double approximation_error;
    PiecewiseSquareLoss model;
    int thj;
    FpopInference(double p, double a, PiecewiseSquareLoss m, int t);
};

FpopInference fpop_analytic_inference_recycle(PiecewiseSquareLosses * cost_model_fwd,
                                       PiecewiseSquareLosses * cost_model_rev,
                                       double * data_vec, int data_count,
                                       double penalty,
                                       int thj,
                                       int window_size,
                                       double sig);

