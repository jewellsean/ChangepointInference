#include "funPieceList.h"
#include <stdexcept>

double min_cost_change_at_thj(PiecewiseSquareLoss *, int, int, double *, double, double, double, double);

double opt_cost_segmenting_change_at_thj(PiecewiseSquareLoss *,
                                         PiecewiseSquareLoss *,
                                         int, // changepoint of interest
                                         int, // size of window around thj
                                         double *, // original data
                                         double *, // flipped data
                                         double, // shift to data, (phi - nu * y) / (2 * window_size)
                                         double, // tuning parameter to penalize the number of spikes
                                         double,
                                         double);


int thj_in_model_at_phi(
        PiecewiseSquareLosses *cost_fwd,
        PiecewiseSquareLosses *cost_rev,
        int thj, // changepoint of interest
        int window_size, // size of window around thj
        int data_count, // number of data points
        double *data_vec, // original data
        double *data_vec_rev, // flipped data
        double theta, // shift to data, (phi - nu * y) / (2 * window_size)
        double penalty, // tuning parameter to penalize the number of spikes
        double *opt_cost, // optimal cost of this segmentation
        int verbose
);


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
                            int verbose);


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
);

double construct_vTy(double * sub_data_f, int n_sub_f, double * sub_data_r, int n_sub_r);
double construct_vTy(int thj, int window_size, double * data_vec, int data_count);
double construct_nu2(int thj, int window_size, int data_count);

void check_selective_inference(PiecewiseSquareLoss * analytic_phi,
                               int thj, // changepoint of interest
                               int window_size, // size of window around thj
                               int data_count, // number of data points
                               double *data_vec, // original data
                               double penalty, // tuning parameter to penalize the number of spikes
                               int verbose);


double calc_p_value(PiecewiseSquareLoss * analytic_phi,
                    int thj, // changepoint of interest
                    int window_size, // size of window around thj
                    int data_count, // number of data points
                    double *data_vec, // original data
                    double sig, // noise variance
                    int verbose);