#pragma once
#include <string>

void print_array(double *ary, int data_count);
void print_int_array(int *ary, int data_count);

double * subset_array(double *data_vec, int start, int end);

double * add_constant_array(double *data_vec, double c, int data_count);

double * reverse_data(double *data_vec, int data_count);

double * read_data(const std::string filename, const int data_count);

double pnorm(double x);
double pnorm_log(double x);
double log_sum_exp(double logx, double logy);
double log1mexp(double a);
double log_subtract(double x, double y);

