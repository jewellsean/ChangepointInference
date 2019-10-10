#include <stdio.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <stdexcept>

void print_array(double *ary, int data_count) {
  printf("[ %f \t", ary[0]);
  for(int data_i = 1; data_i < (data_count - 1); data_i++){
    printf("%f \t", ary[data_i]);
    if (data_i % 10 == 0) {
      printf("\n");
    }
  }
  printf("  %f ]\n", ary[data_count -1]);
}

void print_int_array(int *ary, int data_count) {
  printf("[ %d \t", ary[0]);
  for(int data_i = 1; data_i < (data_count - 1); data_i++){
    printf("%d \t", ary[data_i]);
    if (data_i % 10 == 0) {
      printf("\n");
    }
  }
  printf("  %d ]\n", ary[data_count -1]);
}



/*
 * subset an array, ary[start, end)
 */
double * subset_array(double *data_vec, int start, int end) {
  int n_subset = end - start + 1;
  double *data_subset = new double[n_subset];

  for (int i=start; i < end; i++) {
    data_subset[i - start] = data_vec[i];
  }

  return data_subset;
}


/*
 * add constant to array
 */
double * add_constant_array(double *data_vec, double c, int data_count) {
  double *data_added = new double[data_count];
  for (int i=0; i < data_count; i++) {
    data_added[i] = data_vec[i] + c;
  }
  return data_added;
}


/*
 * reverse a data array
 */
double * reverse_data(double *data_vec, int data_count) {
  double *data_vec_rev = new double[data_count];
  for (int i = 0; i < data_count; i++) {
    data_vec_rev[i] = data_vec[data_count - 1 - i];
  }
  return data_vec_rev;
}


/*
 * read data from a file
 * output as a double * array
 */
double * read_data(const std::string filename, const int data_count) {
  const int max_line_size = 100;
  char in[max_line_size];

  std::ifstream fsnumbers;
  fsnumbers.open(filename);

  double * data_vec = new double[data_count];
  for (int i = 0; i < data_count; i++) {
    fsnumbers >> in;
    data_vec[i] = std::stod(in);
  }
  return data_vec;
}


/*
 * Double precision standard normal cdf
 * https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c/23119456#23119456
 * returns on log scale
 *
 * possibly consider https://github.com/SurajGupta/r-source/blob/master/src/nmath/pnorm.c
 *
 */
double pnorm_log(double x) {
  static const double RT2PI = sqrt(4.0*acos(0.0));

  static const double SPLIT = 7.07106781186547;

  static const double N0 = 220.206867912376;
  static const double N1 = 221.213596169931;
  static const double N2 = 112.079291497871;
  static const double N3 = 33.912866078383;
  static const double N4 = 6.37396220353165;
  static const double N5 = 0.700383064443688;
  static const double N6 = 3.52624965998911e-02;
  static const double M0 = 440.413735824752;
  static const double M1 = 793.826512519948;
  static const double M2 = 637.333633378831;
  static const double M3 = 296.564248779674;
  static const double M4 = 86.7807322029461;
  static const double M5 = 16.064177579207;
  static const double M6 = 1.75566716318264;
  static const double M7 = 8.83883476483184e-02;

  const double z = fabs(x);
  double c = 0.0;

  if(z<=50.0)
  {
    const double e = exp(-z*z/2.0);
    if(z<SPLIT)
    {
      const double n = (((((N6*z + N5)*z + N4)*z + N3)*z + N2)*z + N1)*z + N0;
      const double d = ((((((M7*z + M6)*z + M5)*z + M4)*z + M3)*z + M2)*z + M1)*z + M0;
      c = e*n/d;
    }
    else
    {
      const double f = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
      c = e/(RT2PI*f);
    }
  }
  return x<=0.0 ? log(c) : log1p(-c);
}

double log_sum_exp(double logx, double logy) {

//  printf("\t\t\t\t **** log_sum_exp **** \n");
//  printf("inputs logx = %f, logy = %f \n", logx, logy);
  double a;
  if (logx > logy) {
    a = logx;
  } else {
    a = logy;
  }
  if (fabs(a) == INFINITY) {
    a = 0;
  }

//  printf("maximium value = %f \n", a);

  double out = exp(logx - a);
//  printf("out 1 = %f\n", out);
  out += exp(logy - a);
//  printf("out 2 = %f\n", out);
//  printf("final out = %f\n", log(out) + a);
//  printf("\t\t\t\t **** END OF log_sum_exp **** \n");
  return log(out) + a;
}


// https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
double log1mexp(double a) {
        if (a >= 0 && a <= log(2)) {
          return(log(-expm1(-a)));
        } else if (a > log(2)) {
          return(log1p(-exp(-a)));
        } else {
          std::range_error("log1mexp:: input a must be positive");
        }
}


double log_subtract(double x, double y) {
if (x < y) {
  std::range_error("log_subtract:: cannot take log of (-)ve number");
}
  return(x + log1mexp(fabs(y - x)));
}
