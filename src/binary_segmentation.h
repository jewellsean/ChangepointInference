#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#define ABS(x) ((x)<0 ? -(x) : (x))

#pragma once

typedef std::vector<double> VecDouble;
typedef std::vector<int> VecInt;

class BinarySegmentation {
public:
    VecInt ordered_changepoints;
    VecInt changepoint_signs;
    BinarySegmentation(const VecInt & o, const VecInt & s);
    BinarySegmentation();
    void print();
};


BinarySegmentation k_step_bs(const VecDouble & Sy, const int n_steps);
BinarySegmentation k_step_bs_wrap(const VecDouble & y, const int n_steps);


class CusumInterval {
public:
    double max_csum;
    int argmax_csum;
    int sign_csum;
    CusumInterval(double m, int a, int s);
    CusumInterval();
};

double calculate_csum(const VecDouble & Sy, int tau, int s, int e);