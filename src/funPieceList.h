/* -*- compile-command: "R CMD INSTALL .." -*- */
#pragma once
#include <list>
#include <vector>
#include <stdexcept>
#define NEWTON_EPSILON 1e-30
#define DIFF_EPS 1e-8
#define PREV_NOT_SET (-3)
#define MACHINE_MIN -1e10
#define MACHINE_MAX 1e10
#define TAIL_THRESHOLD 30.0
#define ABS(x) ((x)<0 ? -(x) : (x))

class SquareLossPiece {
public:
  // coefs for square loss
  // Square * u ^ 2 + Linear * u + Constant * 1
  double Square;
  double Linear;
  double Constant;
  double min_mean;
  double max_mean;
  int data_i;
  double prev_mean;
  SquareLossPiece();
  SquareLossPiece
    (double a, double b, double c, double m, double M, int i, double);
  double argmin();
  double argmin_mean();
  void print();
  double get_smaller_root(double);
  double get_larger_root(double);
  bool has_two_roots(double);
  double getCost(double mean);
  double SquareLoss(double);
};


typedef std::list<SquareLossPiece> SquareLossPieceList;

class PiecewiseSquareLoss {
public:
  SquareLossPieceList piece_list;
  void set_to_unconstrained_min_of(PiecewiseSquareLoss *, int);
  void set_to_min_env_of
    (PiecewiseSquareLoss *, PiecewiseSquareLoss *, int);
  void set_to_addition_of
            (PiecewiseSquareLoss *, PiecewiseSquareLoss *, int);
  int check_min_of(PiecewiseSquareLoss *, PiecewiseSquareLoss *);
  void push_min_pieces(PiecewiseSquareLoss *, PiecewiseSquareLoss *, SquareLossPieceList::iterator, SquareLossPieceList::iterator, int);
  void push_piece(SquareLossPieceList::iterator, double, double);
  void add(double Square, double Linear, double Constant);
  void print();
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, double *prev_mean);
  int findSegEnd(double mean);
  double findCost(double mean);
  double getMinCost();
  int getMinIntervalInd();
  void Minimize
    (double *best_cost,
     double *best_mean,
     int *data_i,
     double *prev_mean);
};

double findCostDiff(PiecewiseSquareLoss *f, PiecewiseSquareLoss *g, double mean);
bool sameFunsSquare(SquareLossPieceList::iterator, SquareLossPieceList::iterator);

typedef std::vector<PiecewiseSquareLoss> PiecewiseSquareLosses;
double MidMean(double first, double second);

class BiSquareLossPiece {
public:
    // coefs for bi-square loss in u and p
    // SquareU * u^2 + LinearU * u + LinearUP * u * p + SquareP * p^2 + LinearP * p + Constant
    double SquareU;
    double LinearU;
    double LinearUP;
    double SquareP;
    double LinearP;
    double Constant;
    double min_u;
    double max_u;
    double min_p;
    double max_p;
    BiSquareLossPiece();
    BiSquareLossPiece
            (double, double, double, double, double, double, double, double, double, double);
    void print();
    PiecewiseSquareLoss min_over_u();
};



typedef std::list<BiSquareLossPiece> BiSquareLossPieceList;

class PiecewiseBiSquareLoss {
public:
    BiSquareLossPieceList piece_list;
    void set_to_addition_of
            (PiecewiseBiSquareLoss *, PiecewiseBiSquareLoss *, int);
    void set_to_pw_u(PiecewiseSquareLoss *f);
    void set_to_pw_p(PiecewiseSquareLoss *f);
    void print_quad_p();
    PiecewiseSquareLoss get_univariate_p();
    void add(double, double, double, double, double, double);
    void print();
    PiecewiseSquareLoss min_over_u();
};

typedef std::list<PiecewiseBiSquareLoss> PiecewiseBiSquareLossList;

class PiecewiseBiSquareLosses {
public:
    PiecewiseBiSquareLossList piece_list;
    void add(double, double, double, double, double, double);
    void collect(PiecewiseBiSquareLoss *);
    PiecewiseBiSquareLoss min_u();
    void print();
    void set_to_addition_of(PiecewiseBiSquareLosses *, PiecewiseBiSquareLosses *, int);
};