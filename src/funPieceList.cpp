/* -*- compile-command: "R CMD INSTALL .." -*- */
/* Functions extended from https://github.com/jewellsean/FastLZeroSpikeInference and
 * https://github.com/tdhock/PeakSegOptimal
 * */
#include "funPieceList.h"
#include <list>
#include <math.h>
#include <set>
#include <stdio.h>

#define SGN(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

#define SET_DIFF(x) ((ABS(x) > DIFF_EPS) ? x : 0)

SquareLossPiece::SquareLossPiece
  (double a, double b, double c, double m, double M, int i, double prev){
  Square = a;
  Linear = b;
  Constant = c;
  min_mean = m;
  max_mean = M;
  data_i = i;
  prev_mean = prev;
}

SquareLossPiece::SquareLossPiece(){
}

bool SquareLossPiece::has_two_roots(double equals){
  // are there two solutions to the equation 
  // Square * u ^ 2 + 
  // Linear * u + Constant = equals ? 
  double delta = Linear * Linear - 4 * Square * (Constant - equals);
  if (delta > 0) {
    return true;
  } else {
    return false;
  }
}

double SquareLossPiece::SquareLoss(double mean){
  return Square * mean * mean + Linear * mean + Constant;
}

double log_discriminant(double A, double B, double C, int sa, int sc) {
// A = log(|a|), B = log(|b|), C = log(|c|)
// sa = sign(a), sb = sign(b), sc = sign(c)
  double max_val = std::max(std::max(A, B), C);
  if (sa * sc > 0) {
    return max_val + log(exp(2 * B - max_val) - exp(log(4) + A + C - max_val));
  } else {
    return max_val + log(exp(2 * B - max_val) + exp(log(4) + A + C - max_val));
  }
}


/*
  diff_of_products() computes a*b-c*d with a maximum error < 1.5 ulp

  Claude-Pierre Jeannerod, Nicolas Louvet, and Jean-Michel Muller,
  "Further Analysis of Kahan's Algorithm for the Accurate Computation
  of 2x2 Determinants". Mathematics of Computation, Vol. 82, No. 284,
  Oct. 2013, pp. 2245-2264

  https://stackoverflow.com/questions/48979861/numerically-stable-method-for-solving-quadratic-equations/50065711

  Also see: https://www.boost.org/doc/libs/1_70_0/boost/math/tools/roots.hpp

*/
double diff_of_products (double a, double b, double c, double d)
{
  double w = d * c;
  double e = fma (-d, c, w);
  double f = fma (a, b, -w);
  return f + e;
}



double SquareLossPiece::get_larger_root(double equals){
  double q = -0.5 * (Linear + copysign (sqrt (diff_of_products (Linear, Linear, 4.0*Square, Constant - equals)), Linear));
  double r1 = q / Square;
  double r2 = (Constant - equals) / q;
  if (r1 > r2) {
    return r1;
  } else {
    return r2;
  }
}

double SquareLossPiece::get_smaller_root(double equals){
  double q = -0.5 * (Linear + copysign (sqrt (diff_of_products (Linear, Linear, 4.0*Square, Constant - equals)), Linear));
  double r1 = q / Square;
  double r2 = (Constant - equals) / q;
  if (r1 < r2) {
    return r1;
  } else {
    return r2;
  }
}

double SquareLossPiece::argmin_mean(){
  // f(u) = Square * u ^ 2 + Linear * u + Constant,
  // f'(u)= 2 * Square * u + Linear = 0 means
  // u = -Linear / (2 * Square)
  
  if (Square > 0 || Square < 0) {
    return - Linear / (2 * Square);  
  } else if((Square == 0) & (Linear > 0)) {
    return min_mean;
  } else if((Square == 0) & (Linear < 0)) {
    return max_mean; 
  } else if((Square == 0) & (Linear == 0)) {
   return min_mean; 
  }
  throw 1;
}

double SquareLossPiece::argmin(){
  return argmin_mean();
}

double SquareLossPiece::getCost(double mean){
    if (mean < INFINITY && mean > -INFINITY) {
      return Square * mean * mean + Linear * mean + Constant;
    } else if (mean == INFINITY) {
      if (Square > 0) {
        return INFINITY;
      } else if (Square < 0) {
        return -INFINITY;
      } else if (Square == 0) {
        if (Linear > 0) {
          return INFINITY;
        } else if (Linear < 0) {
            return -INFINITY;
          } else if (Linear == 0) {
            return Constant;
          }
        }
  } else if (mean == -INFINITY) {
      if (Square > 0) {
        return INFINITY;
      } else if (Square < 0) {
        return -INFINITY;
      } else if (Square == 0) {
        if (Linear > 0) {
          return -INFINITY;
        } else if (Linear < 0) {
          return INFINITY;
        } else if (Linear == 0) {
          return Constant;
        }
      }
    }
  throw std::runtime_error("cannot determine cost of mean");
}


void PiecewiseSquareLoss::set_to_unconstrained_min_of
  (PiecewiseSquareLoss *input, int verbose) {
  piece_list.clear();
  SquareLossPieceList::iterator it = input->piece_list.begin();

  double prev_min_cost = INFINITY;
  double left_most = INFINITY;
  double right_most = -INFINITY;
  double right_mean, left_mean, prev_best_mean, mu_cost, mu;
  
  while(it != input->piece_list.end()){
    if (verbose) {
      printf("start new iter of set to unconstrained min of--------------\n");
      printf("Searching for min cost in \n");
      printf("%10s %10s %15s %15s %15s %15s %s\n",
             "Square", "Linear", "Constant",
             "min_mean", "max_mean",
             "prev_mean", "data_i");
      it -> print();

    }
      
    // boundary costs
    left_mean = it -> min_mean;
    
    if (left_mean < left_most) {
      left_most = left_mean;
    }
    
    right_mean = it -> max_mean;
    
    if (right_mean > right_most) {
      right_most = right_mean;
    }
    
    // determine argmin of this segment
    mu = it->argmin();
    // argmin occurs in interval
    if (mu >= left_mean && mu <= right_mean) {
      mu_cost = it->getCost(mu);
    } else {
      double left_cost = it -> getCost(left_mean);
      double right_cost = it -> getCost(right_mean);
      
      if (right_cost < left_cost) {
        mu_cost = right_cost;
        mu = right_mean;
      } else {
        mu_cost = left_cost;
        mu = left_mean;
      }
    }
    
    if (mu_cost < prev_min_cost) {
      prev_min_cost = mu_cost;
      prev_best_mean = mu;
    }
    it++;
  } // end loop over components
  
    piece_list.emplace_back
    (0, 0, prev_min_cost,
     left_most, right_most, PREV_NOT_SET,
     prev_best_mean);
  
  if (verbose) {
    printf("interval [%f, %f]\n", left_most, right_most);
    printf("Minimum cost %f \n", prev_min_cost);
    printf("------------------------------------------\n");
  }
}

void PiecewiseSquareLoss::add(double Square, double Linear, double Constant){
  SquareLossPieceList::iterator it;
    for(it=piece_list.begin(); it != piece_list.end(); it++){
      it->Square += Square;
      it->Linear += Linear;
      it->Constant += Constant;
    } 
  }


void PiecewiseSquareLoss::set_prev_seg_end(int prev_seg_end){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->data_i = prev_seg_end;
  }
}

void PiecewiseSquareLoss::findMean
  (double mean, int *seg_end, double *prev_mean){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    //	  printf("looking for mean %f in interval [%f, %f]\n", mean, it -> min_mean, it -> max_mean);
    if(it->min_mean <= mean && mean <= it->max_mean){
      *seg_end = it->data_i;
      *prev_mean = it->prev_mean;
      //      printf("peel off mean %f \n", &(it -> prev_mean));
      return;
    }
  }
}

int PiecewiseSquareLoss::findSegEnd(double mean){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    //	  printf("looking for mean %f in interval [%f, %f]\n", mean, it -> min_mean, it -> max_mean);
    if(it->min_mean <= mean && mean <= it->max_mean){
      return it->data_i;
    }
  }
}


double PiecewiseSquareLoss::findCost(double mean){
  SquareLossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_mean <= mean && mean <= it->max_mean){
      return it->getCost(mean);
    }
  }
  throw std::range_error("mean not contained in piecewise list");
}

double findCostDiff(PiecewiseSquareLoss *f, PiecewiseSquareLoss *g, double mean){
  SquareLossPieceList::iterator it1, it2;
  it1 = f->piece_list.begin();
  it2 = g->piece_list.begin();

  while (it1 != f->piece_list.end() && !(it1 -> min_mean <= mean && it1 -> max_mean > mean)) it1++;
  while (it2 != g->piece_list.end() && !(it2 -> min_mean <= mean && it2 -> max_mean > mean)) it2++;

  SquareLossPiece diff_piece
          (it1->Square - it2->Square,
           it1->Linear - it2->Linear,
           it1->Constant - it2->Constant,
           std::max(it1 -> min_mean, it2 -> min_mean),
           std::min(it1 -> max_mean, it2 -> max_mean),
           -5, false);
  return(diff_piece.getCost(mean));

}

void PiecewiseSquareLoss::print(){
  SquareLossPieceList::iterator it;
  printf("%10s %10s %15s %15s %15s %15s %s\n",
         "Square", "Linear", "Constant",
         "min_mean", "max_mean",
         "prev_mean", "data_i");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}


void SquareLossPiece::print(){
  printf("%.15f %.15f %.15f %.5f %.5f %.5f %d\n",
         Square, Linear, Constant,
         min_mean, max_mean,
         prev_mean, data_i);
}

void PiecewiseSquareLoss::Minimize
  (double *best_cost,
   double *best_mean,
   int *data_i,
   double *prev_mean){
  double candidate_cost, candidate_mean;
  int verbose=false;
  SquareLossPieceList::iterator it;
  *best_cost = INFINITY;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_mean = it->argmin();
    if(candidate_mean < it->min_mean){
      candidate_mean = it->min_mean;
    }else if(it->max_mean < candidate_mean){
      candidate_mean = it->max_mean;
    }
    candidate_cost = it->getCost(candidate_mean);
    if(candidate_cost < *best_cost){
      *best_cost = candidate_cost;
      *best_mean = candidate_mean;
      *data_i = it->data_i;
      *prev_mean = it->prev_mean;
    }
  }
}

double PiecewiseSquareLoss::getMinCost() {
  double candidate_cost, candidate_mean;
  int verbose=false;
  SquareLossPieceList::iterator it;
  double best_cost = INFINITY;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_mean = it->argmin();
    if(candidate_mean < it->min_mean){
      candidate_mean = it->min_mean;
    }else if(it->max_mean < candidate_mean){
      candidate_mean = it->max_mean;
    }
    candidate_cost = it->getCost(candidate_mean);
    if(candidate_cost < best_cost){
      best_cost = candidate_cost;
    }
  }
  return best_cost;
}


int PiecewiseSquareLoss::getMinIntervalInd() {
  double candidate_cost, candidate_mean;
  int verbose=false;
  SquareLossPieceList::iterator it;
  double best_cost = INFINITY;
  int data_best;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    candidate_mean = it->argmin();
    if(candidate_mean < it->min_mean){
      candidate_mean = it->min_mean;
    }else if(it->max_mean < candidate_mean){
      candidate_mean = it->max_mean;
    }
    candidate_cost = it->getCost(candidate_mean);
    if(candidate_cost < best_cost){
      best_cost = candidate_cost;
      data_best = it->data_i;
    }
  }
  return data_best;
}

double MidMean(double first, double second) {
  if (first == -INFINITY && second == INFINITY) {
    return 0.0;
  } else if (first == -INFINITY && second != INFINITY) {
    return second - 1;
  } else if (first != -INFINITY && second == INFINITY) {
    return first + 1;
  } else {
    return ((first + second)  / 2);
  }
}

// check that this function is the minimum on all pieces.
// convection is that intervals are of form (a, b]
int PiecewiseSquareLoss::check_min_of
  (PiecewiseSquareLoss *prev, PiecewiseSquareLoss *model){
  SquareLossPieceList::iterator it;
  double cost_diff;
  for(it = piece_list.begin(); it != piece_list.end(); it++){
    if(it != piece_list.begin()){
      SquareLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
        printf("prev->max_mean != it->min_mean min\n");
        return 3;
      }
    }
    if(it->max_mean -  it->min_mean <= -NEWTON_EPSILON){
      printf("max_mean (=%15.10f) <=min_mean(=%15.10f) min\n", it -> max_mean, it->min_mean);
      return 2;
    }
    double mid_mean = MidMean(it -> min_mean, it -> max_mean);
    if(-INFINITY < mid_mean){
      double cost_min = it->getCost(mid_mean);
      double cost_prev = prev->findCost(mid_mean);

      cost_diff = findCostDiff(prev, this, mid_mean);

      if(cost_diff < -DIFF_EPS){
        printf("cost diff(%f)=%15.20f\n", mid_mean, cost_diff);
        printf("prev(%f)=%f\n", mid_mean, cost_prev);
        prev->print();
        printf("min(%f)=%f\n", mid_mean, cost_min);
        print();
        return 1;
      }
      double cost_model = model->findCost(mid_mean);
      cost_diff = findCostDiff(model, this, mid_mean);
      if(cost_diff < -DIFF_EPS){
        printf("cost diff(%f)=%15.20f\n", mid_mean, cost_diff);
        printf("model(%f)=%f\n", mid_mean, cost_model);
        model->print();
        printf("min(%f)=%f\n", mid_mean, cost_min);
        print();
        return 1;
      }
    }
  }
  for(it = prev->piece_list.begin(); it != prev->piece_list.end(); it++){
    if(it != prev->piece_list.begin()){
      SquareLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
        printf("prev->max_mean != it->min_mean prev\n");
        return 3;
      }
    }
    if(it->max_mean - it->min_mean <= -NEWTON_EPSILON){
      printf("max_mean<=min_mean=%15.10f prev\n", it->min_mean);
      return 2;
    }
    double mid_mean = MidMean(it -> min_mean, it -> max_mean);
    if(-INFINITY < mid_mean){
      double cost_prev = it->getCost(mid_mean);
      double cost_min = findCost(mid_mean);

      cost_diff = findCostDiff(prev, this, mid_mean);

      if(cost_diff < -DIFF_EPS){
        printf("cost diff(%f)=%15.20f\n", mid_mean, cost_diff);
        printf("prev(%f)=%f\n", mid_mean, cost_prev);
        prev->print();
        printf("min(%f)=%f\n", mid_mean, cost_min);
        print();
        return 1;
      }
    }
  }
  for(it = model->piece_list.begin(); it != model->piece_list.end(); it++){
    if(it != model->piece_list.begin()){
      SquareLossPieceList::iterator pit = it;
      pit--;
      if(pit->max_mean != it->min_mean){
        printf("prev->max_mean != it->min_mean model\n");
        return 3;
      }
    }
    if(it->max_mean - it->min_mean <= -NEWTON_EPSILON){
      printf("max_mean=%15.10f<=min_mean=%15.10f model\n", it-> max_mean, it->min_mean);
      return 2;
    }
    double mid_mean = MidMean(it -> min_mean, it -> max_mean);
    if(-INFINITY < mid_mean){
      double cost_model = it->getCost(mid_mean);
      double cost_min = findCost(mid_mean);

      cost_diff = findCostDiff(model, this, mid_mean);

      if(cost_diff < -DIFF_EPS){
        printf("cost diff(%f)=%15.20f\n", mid_mean, cost_diff);
        printf("model(%f)=%f\n", mid_mean, cost_model);
        model->print();
        printf("min(%f)=%f\n", mid_mean, cost_min);
        print();
        return 1;
      }
    }
  }
  return 0;
}

void PiecewiseSquareLoss::set_to_min_env_of
  (PiecewiseSquareLoss *fun1, PiecewiseSquareLoss *fun2, int verbose){
  SquareLossPieceList::iterator
  it1 = fun1->piece_list.begin(),
  it2 = fun2->piece_list.begin();
  if(verbose){
    printf("computing min env of:\n");
    printf("=min-less/more\n");
    fun1->print();
    printf("=cost model\n");
    fun2->print();
  }
  piece_list.clear();

  while(it1 != fun1->piece_list.end() &&
        it2 != fun2->piece_list.end()){
    push_min_pieces(fun1, fun2, it1, it2, verbose);
    if(verbose){
      print();
      printf("------\n");
    }
    double last_max_mean = piece_list.back().max_mean;
    if(it1->max_mean == last_max_mean){
      it1++;
    }
    if(it2->max_mean == last_max_mean){
      it2++;
    }
  }
}

bool sameFunsSquare(SquareLossPieceList::iterator it1,
   SquareLossPieceList::iterator it2){
  return it1->Linear == it2->Linear &&
    it1->Square == it2->Square &&
    ABS(it1->Constant - it2->Constant) < NEWTON_EPSILON;
}

void PiecewiseSquareLoss::push_min_pieces
        (PiecewiseSquareLoss *fun1,
         PiecewiseSquareLoss *fun2,
         SquareLossPieceList::iterator it1,
         SquareLossPieceList::iterator it2,
         int verbose){
  bool same_at_left;
  double last_min_mean;
  SquareLossPieceList::iterator prev2 = it2;
  prev2--;
  SquareLossPieceList::iterator prev1 = it1;
  prev1--;


  if (verbose) {
    printf("it1 = \n");
    it1 -> print();
    printf("it2 = \n");
    it2 -> print();
  }


  if(it1->min_mean < it2->min_mean){
    //it1 function piece starts to the left of it2.
    same_at_left = sameFunsSquare(prev2, it1);
    last_min_mean = it2->min_mean;
  }else{
    //it1 function piece DOES NOT start to the left of it2.
    last_min_mean = it1->min_mean;
    if(it2->min_mean < it1->min_mean){
      //it2 function piece starts to the left of it1.
      same_at_left = sameFunsSquare(prev1, it2);
    }else{
      //it1 and it2 start at the same min_mean value.
      if(it1==fun1->piece_list.begin() &&
         it2==fun2->piece_list.begin()){
        same_at_left = false;
      }else{
        same_at_left = sameFunsSquare(prev1, prev2);
      }
    }
  }
  SquareLossPieceList::iterator next2 = it2;
  next2++;
  SquareLossPieceList::iterator next1 = it1;
  next1++;
  bool same_at_right;
  double first_max_mean;
  if(it1->max_mean < it2->max_mean){
    if(verbose)printf("it2 function piece continues to the right of it1.\n");
    same_at_right = sameFunsSquare(next1, it2);
    first_max_mean = it1->max_mean;
  }else{
    first_max_mean = it2->max_mean;
    if(it2->max_mean < it1->max_mean){
      if(verbose)printf("it2 function piece ends before it1.\n");
      same_at_right = sameFunsSquare(it1, next2);
    }else{
      if(verbose)printf("it2 and it1 end at same max_mean.\n");
      if(next1==fun1->piece_list.end() &&
         next2==fun2->piece_list.end()){
        if(verbose)printf("at the end so they can't be equal after this interval.\n");
        same_at_right = false;
      }else{
        if(verbose){
          printf("comparing next function pieces.\n");
          next1->print();
          next2->print();
        }
        same_at_right = sameFunsSquare(next1, next2);
      }
    }
  }
  if(last_min_mean == first_max_mean){
    // we should probably never get here, but if we do, no need to
    // store this interval.
    if(verbose){
      printf("prev\n");
      fun1->print();
      printf("model\n");
      fun2->print();
      printf("interval size 0!-----------------\n");
    }

    // if ((it1 -> Constant) < (it2 -> Constant)) {
    //   push_piece(it1, last_min_mean, first_max_mean);
    // } else {
    //   push_piece(it2, last_min_mean, first_max_mean);
    // }
    return;
  }
  if(sameFunsSquare(it1, it2)){
    // The functions are exactly equal over the entire interval so we
    // can push either of them.
    push_piece(it1, last_min_mean, first_max_mean);
    if(verbose)printf("exactly equal over entire interval\n");
    return;
  }
  SquareLossPiece diff_piece
          (it1->Square - it2->Square,
           it1->Linear - it2->Linear,
           it1->Constant - it2->Constant,
           last_min_mean, first_max_mean,
           -5, false);
  // Evaluate the middle in the original space, to avoid problems when
  // SWJ: need to be careful here in case either first_max_mean of last_min_mean
  // is -Inf or Inf
  // If either is Inf or -Inf then

  double mid_mean;
  if (first_max_mean != INFINITY && last_min_mean != -INFINITY) {
    mid_mean = (first_max_mean + last_min_mean) / 2;
  } else if (first_max_mean != INFINITY && last_min_mean == -INFINITY) {
    mid_mean = first_max_mean - 1;
  } else if (first_max_mean == INFINITY && last_min_mean != -INFINITY) {
    mid_mean = last_min_mean + 1;
  } else  {
    mid_mean = 0;
  }

  double cost_diff_mid = diff_piece.getCost(mid_mean);
  // Easy case of equality on both left and right.
  if(same_at_left && same_at_right){
    if(verbose)printf("Same on both the left and the right\n");
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    return;
  }

  // add back the easy degenerate cases that do not require root finding

  // Easy degenerate cases that do not require root finding.
  if(diff_piece.Square == 0){
    // g(x) = Linear * x + Constant = 0,
    // x = - Constant / Linear
    if(diff_piece.Linear == 0){
      // They are offset by a Constant.
      if(diff_piece.Constant < 0){
        push_piece(it1, last_min_mean, first_max_mean);
      }else{
        push_piece(it2, last_min_mean, first_max_mean);
      }
      if(verbose)printf("offset by a constant=%e\n", diff_piece.Constant);
      return;
    }
    if(diff_piece.Constant == 0){
      // The only difference is the Linear coef.
      if(cost_diff_mid < 0){
        push_piece(it1, last_min_mean, first_max_mean);
      }else{
        push_piece(it2, last_min_mean, first_max_mean);
      }
      if(verbose)printf("only diff is linear coef\n");
      return;
    }
    double mean_at_equal_cost = -diff_piece.Constant/diff_piece.Linear;
    if(last_min_mean < mean_at_equal_cost &&
       mean_at_equal_cost < first_max_mean){
      // the root is in the interval, so we need to add two intervals.
      if(diff_piece.Linear > 0){
        push_piece(it1, last_min_mean, mean_at_equal_cost);
        push_piece(it2, mean_at_equal_cost, first_max_mean);
      }else{
        push_piece(it2, last_min_mean, mean_at_equal_cost);
        push_piece(it1, mean_at_equal_cost, first_max_mean);
      }
      if(verbose)printf("Square zero with one root in interval\n");
      return;
    }
    // the root is outside the interval, so one is completely above
    // the other over this entire interval.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    if(verbose)printf("Square zero with no roots in interval\n");
    return;
  }


  double cost_diff_left = diff_piece.getCost(last_min_mean);
  double cost_diff_right = diff_piece.getCost(first_max_mean);
  bool two_roots = diff_piece.has_two_roots(0.0);

//  printf("difference in left and mid points %.50f\n", last_min_mean - mid_mean);
//  printf("cost at left %0.50f\n", diff_piece.getCost(last_min_mean));
//  printf("cost at mid %0.50f\n", diff_piece.getCost(mid_mean));
//  printf("cost diff between left and mid = %.50f\n", diff_piece.getCost(last_min_mean) - diff_piece.getCost(mid_mean));
//
//  printf("smaller root at 0 cost %f\n", diff_piece.get_smaller_root(0.0));
//  printf("larger root at 0 cost %f\n", diff_piece.get_larger_root(0.0));
//  printf("discriminant %0.50f\n", diff_piece.Linear * diff_piece.Linear - 4 * diff_piece.Square * diff_piece.Constant);
//
//
//  printf("----\n");
//  printf("manual cost calc @ left %0.50f\n", diff_piece.Square * last_min_mean * last_min_mean + diff_piece.Linear * last_min_mean + diff_piece.Constant);
//  printf("manual cost calc @ mid point %0.50f\n", diff_piece.Square * mid_mean * mid_mean + diff_piece.Linear * mid_mean + diff_piece.Constant);
//  printf("mid point in defined interval %d\n", diff_piece.min_mean <= mid_mean && mid_mean <= diff_piece.max_mean);

  double smaller_mean, larger_mean;
  if(two_roots){
    smaller_mean = diff_piece.get_smaller_root(0.0);
    larger_mean = diff_piece.get_larger_root(0.0);
  }
  if(same_at_right){
    // they are equal on the right, but we don't know if there is
    // another crossing point somewhere to the left.
    if(two_roots){
      // there could be a crossing point to the left.
      double mean_at_crossing = smaller_mean;
      double mean_between_zeros = (mean_at_crossing + first_max_mean)/2;
      double cost_between_zeros = diff_piece.getCost(mean_between_zeros);
      double mean_at_optimum = diff_piece.argmin();
      if(verbose){
        printf("cost_diff(left:%e)=%e\n", last_min_mean, cost_diff_left);
        printf("cost_diff(cross:%e)=%e\n", mean_at_crossing, diff_piece.getCost(mean_at_crossing));
        printf("cost_diff(between:%e)=%e\n", mean_between_zeros, cost_between_zeros);
        printf("cost_diff(optimum:%e)=%e\n", mean_at_optimum, diff_piece.getCost(mean_at_optimum));
        printf("cost_diff(right:%e)=%e\n", first_max_mean, cost_diff_right);
      }
      if(last_min_mean < mean_at_crossing &&
         mean_at_crossing < mean_at_optimum &&
         mean_at_optimum < first_max_mean){
        //the cross point is in the interval.
        if(cost_diff_left < 0){
          push_piece(it1, last_min_mean, mean_at_crossing);
          push_piece(it2, mean_at_crossing, first_max_mean);
        }else{
          push_piece(it2, last_min_mean, mean_at_crossing);
          push_piece(it1, mean_at_crossing, first_max_mean);
        }
        if(verbose)printf("equal on the right with one crossing in interval\n");
        return;
      }
    }//if(two_roots
    // Test the cost at the midpoint, since the cost may be equal on
    // both the left and the right.
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    if(verbose)printf("equal on the right with no crossing in interval\n");
    return;
  }
  if(same_at_left){
    // equal on the left.
    if(two_roots){
      // There could be a crossing point to the right.
      double mean_at_crossing = larger_mean;
      double mean_at_optimum = diff_piece.argmin();
      if(verbose)printf("mean_at_crossing=%f\n", mean_at_crossing);
      if(verbose)printf("last_min_mean=%f\n", last_min_mean);
      if(verbose)printf("mean_at_optimum=%f\n", mean_at_optimum);
      if(verbose)printf("first_max_mean=%f\n", first_max_mean);
      if(last_min_mean < mean_at_optimum &&
         mean_at_optimum < mean_at_crossing &&
         mean_at_crossing < first_max_mean){
        // the crossing point is in this interval.
        if(cost_diff_right < 0){
          push_piece(it2, last_min_mean, mean_at_crossing);
          push_piece(it1, mean_at_crossing, first_max_mean);
        }else{
          push_piece(it1, last_min_mean, mean_at_crossing);
          push_piece(it2, mean_at_crossing, first_max_mean);
        }
        if(verbose)printf("equal on the left with crossing in interval\n");
        return;
      }
    }//if(there may be crossing
    if(cost_diff_mid < 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
    // printf("same at left value %d \n", same_at_left);
    if (verbose) printf("cost diff left, mid, right (%f, %f, %f) \n", cost_diff_left, cost_diff_mid, cost_diff_right);
    if(verbose)printf("equal on the left with no crossing in interval\n");
    return;
  }
  // The only remaining case is that the curves are equal neither on
  // the left nor on the right of the interval. However they may be
  // equal inside the interval, so let's check for that.
  double first_mean = INFINITY, second_mean = INFINITY;
  if(two_roots){
    bool larger_inside =
            last_min_mean < larger_mean && larger_mean < first_max_mean;
    if(verbose)printf("smaller_mean=%f \nlarger_mean=%f \n",
                      smaller_mean,
                      larger_mean);
    bool smaller_inside =
            last_min_mean < smaller_mean &&
            // 0 < smaller_mean && // breaks for negative data
            smaller_mean < first_max_mean;
    if(larger_inside){
      if(smaller_inside && smaller_mean < larger_mean){
        // both are in the interval.
        first_mean = smaller_mean;
        second_mean = larger_mean;
        if(verbose){
          diff_piece.print();
          printf("%f and %f in [%f,%f]\n",
                 smaller_mean, larger_mean,
                 last_min_mean, first_max_mean);
        }
      }else{
        // smaller mean is not in the interval, but the larger is.
        first_mean = larger_mean;
        if(verbose){
          printf("%f in [%f,%f]\n",
                 first_mean,
                 last_min_mean, first_max_mean);
        }
      }
    }else{
      // larger mean is not in the interval
      if(smaller_inside){
        // smaller mean is in the interval, but not the larger.
        first_mean = smaller_mean;
        if(verbose){
          printf("%f in [%f,%f]\n",
                 first_mean,
                 last_min_mean, first_max_mean);
        }
      }
    }
  }//if(two_roots
  if(second_mean != INFINITY){
    // two crossing points.
    double squareDiff = diff_piece.Square;
    if(squareDiff < 0){
      push_piece(it1, last_min_mean, first_mean);
      push_piece(it2, first_mean, second_mean);
      push_piece(it1, second_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_mean);
      push_piece(it1, first_mean, second_mean);
      push_piece(it2, second_mean, first_max_mean);
    }
    if(verbose)printf("not equal on the sides, 2 crossing points\n");
  }else if(first_mean != INFINITY){
    // "one" crossing point.
    // need to be careful here, too
    // last_min_mean could be -Inf
    double before_mean;
    if (last_min_mean != -INFINITY) {
      before_mean = (last_min_mean + first_mean) / 2;
    } else {
      before_mean = first_mean - 1;
    }

    double cost_diff_before = diff_piece.getCost(before_mean);
    if(verbose){
      printf("cost_diff_before(%.55f)=%f\n", before_mean, cost_diff_before);
    }

    // SWJ: need to be careful as first_max_mean could be inf
    double after_mean;

    if (first_max_mean != INFINITY) {
      after_mean = (first_max_mean + first_mean)/2;
    } else {
      after_mean = first_mean + 1; // should be able to check anywhere in the interval (first_mean, inf)?
    }


    double cost_diff_after = diff_piece.getCost(after_mean);
    if(verbose)printf("cost_diff_after(%.55f)=%f\n", after_mean, cost_diff_after);
    if(cost_diff_before < 0){
      if(cost_diff_after < 0){
        // f1-f2<0 meaning f1<f2 on the entire interval, so just push it1.
        push_piece(it1, last_min_mean, first_max_mean);
      }else{
        push_piece(it1, last_min_mean, first_mean);
        push_piece(it2, first_mean, first_max_mean);
      }
    }else{//f1(before)-f2(before)>=0 meaning f1(before)>=f2(before)
      if(cost_diff_after < 0){
        //f1(after)-f2(after)<0 meaning f1(after)<f2(after)
        push_piece(it2, last_min_mean, first_mean);
        push_piece(it1, first_mean, first_max_mean);
      }else{
        //f1(after)-f2(after)>=0 meaning f1(after)>=f2(after)
        push_piece(it2, last_min_mean, first_max_mean);
      }
    }
    if(verbose)printf("not equal on the sides, 1 crossing point\n");
  }else{
    // "zero" crossing points.
    if(verbose){
      printf("not equal on the sides, zero crossing points\n");
      printf("cost_diff left=%e mid=%e right=%e\n",
             cost_diff_left, cost_diff_mid, cost_diff_right);
    }
    double cost_diff;
    if(first_max_mean == INFINITY){
      double minimizer = diff_piece.argmin();
      bool min_inside =
              last_min_mean < minimizer && minimizer < first_max_mean;
      if (min_inside) {
        cost_diff = diff_piece.getCost(minimizer);
      } else {
        cost_diff = diff_piece.getCost(last_min_mean + 1);
      }

    }else{
      if(ABS(cost_diff_mid) < NEWTON_EPSILON){
        cost_diff = cost_diff_right;
      }else{
        cost_diff = cost_diff_mid;
      }
    }
    if(cost_diff <= 0){
      push_piece(it1, last_min_mean, first_max_mean);
    }else{
      push_piece(it2, last_min_mean, first_max_mean);
    }
  }
}

void PiecewiseSquareLoss::push_piece(SquareLossPieceList::iterator it, double min_mean, double max_mean){
  if(max_mean <= min_mean){
    return;
  }
  SquareLossPieceList::iterator last=piece_list.end();
  --last;
  if(piece_list.size() &&
  sameFunsSquare(last, it) &&
  it->prev_mean == last->prev_mean &&
  it->data_i == last->data_i){
    //it is the same function as last, so just make last extend
    //farther to the right.
    last->max_mean = max_mean;
  }else{
    //it is a different function than last, so push it on to the end
    //of the list.
    piece_list.emplace_back(it->Square, it->Linear, it->Constant,
     min_mean, max_mean,
     it->data_i, it->prev_mean);
  }
}

void PiecewiseSquareLoss::set_to_addition_of(PiecewiseSquareLoss *fun1, PiecewiseSquareLoss *fun2, int verbose) {
  SquareLossPieceList::iterator
          it1 = fun1->piece_list.begin(),
          it2 = fun2->piece_list.begin();
  if(verbose){
    printf("computing addition of:\n");
    printf("=function 1\n");
    fun1->print();
    printf("=function 2\n");
    fun2->print();
  }

  piece_list.clear();

  double last_max_mean;
  double max_mean;
  if (it1 -> min_mean > it2 -> min_mean) {
    last_max_mean = it1 -> min_mean;
  } else {
    last_max_mean = it2 -> min_mean;
  }


  while(it1 != fun1->piece_list.end() &&
        it2 != fun2->piece_list.end()){

    if (it1 -> max_mean < it2 -> max_mean) {
      max_mean = it1 -> max_mean;
    } else {
      max_mean = it2 -> max_mean;
    }


    if(verbose) {
      printf("adding from [min_mean, max_mean] = [%f, %f]\n", last_max_mean, max_mean);
      it1->print();
      it2->print();
    }


    piece_list.emplace_back(it1->Square + it2->Square, it1->Linear + it2->Linear, it1->Constant + it2->Constant,
                            last_max_mean, max_mean,
                            NULL, NULL);

    if(it1->max_mean == max_mean){
      it1++;
    }
    if(it2->max_mean == max_mean){
      it2++;
    }
    last_max_mean = max_mean;
  } // end while

}



BiSquareLossPiece::BiSquareLossPiece(double squareU,
        double linearU,
        double linearUP,
        double squareP,
        double linearP,
        double constant,
        double m_u,
        double M_u,
        double m_p,
        double M_p){

SquareU = squareU;
LinearU = linearU;
LinearUP = linearUP;
SquareP = squareP;
LinearP = linearP;
Constant = constant;
min_u = m_u;
max_u = M_u;
min_p = m_p;
max_p = M_p;
}

BiSquareLossPiece::BiSquareLossPiece(){
}

PiecewiseSquareLoss BiSquareLossPiece::min_over_u() {
  /*
   * return piecewise square loss over p
   * each piece produces three intervals in phi
   *
   *
   */
  PiecewiseSquareLoss out;
  out.piece_list.clear();

  if (SquareU == 0 && LinearU == 0 && Constant == 0) { // degenerate case
    out.piece_list.emplace_back(0, 0, 0, min_p, max_p, PREV_NOT_SET, -3);
    return out;
  }

  if (SquareU <= 0) {
    printf("Error with this function \n");
    print();
    throw std::range_error("SquareU is <= 0.");
  }

  if (LinearUP == 0) {
    double constant;
    if (- LinearU / (2 * SquareU) < min_u) {
      constant = min_u * min_u * SquareU + min_u * LinearU + Constant;
    } else if (- LinearU / (2 * SquareU) > max_u) {
      constant = max_u * max_u * SquareU + max_u * LinearU + Constant;
    } else {
      constant =  - LinearU * LinearU / (4 * SquareU) + Constant;
    }

    if (min_p > MACHINE_MIN) {
      out.piece_list.emplace_back(0, 0, INFINITY, MACHINE_MIN, min_p, PREV_NOT_SET, -3);
    }
    out.piece_list.emplace_back(SquareP, LinearP, constant, min_p, max_p, PREV_NOT_SET, -3);
    if (max_p < MACHINE_MAX) {
      out.piece_list.emplace_back(0, 0, INFINITY, max_p, MACHINE_MAX, PREV_NOT_SET, -3);
    }


  } else {
    double end_1, end_2, region_1, region_2;
    if (LinearUP > 0) {
      end_1 = -(2 * max_u * SquareU + LinearU) / LinearUP;
      end_2 = -(2 * min_u * SquareU + LinearU) / LinearUP;
      region_1 = max_u;
      region_2 = min_u;
    } else {
      end_2 = -(2 * max_u * SquareU + LinearU) / LinearUP;
      end_1 = -(2 * min_u * SquareU + LinearU) / LinearUP;
      region_2 = max_u;
      region_1 = min_u;
    }

      if (min_p > MACHINE_MIN) {
        out.piece_list.emplace_back(0, 0, INFINITY, MACHINE_MIN, min_p, PREV_NOT_SET, -3);
      }
      if (min_p < end_1 && max_p >= end_1) {
        out.piece_list.emplace_back(SquareP, region_1 * LinearUP + LinearP,
                                    region_1 * region_1 * SquareU + region_1 * LinearU + Constant, min_p, end_1, -3, PREV_NOT_SET);
      }
      if (max_p < end_1) {
        out.piece_list.emplace_back(0, 0, INFINITY, min_p, max_p, PREV_NOT_SET, -3);
      }
      double start_interval = (min_p > end_1 ? min_p : end_1);
      double end_interval = (max_p > end_2 ? end_2 : max_p);
      if ((end_interval -  start_interval) > 0) {
        out.piece_list.emplace_back(SquareP - LinearUP * LinearUP / (4 * SquareU),
                                    - LinearU * LinearUP / (2 * SquareU) + LinearP,
                                    Constant - LinearU * LinearU / (4 * SquareU), start_interval, end_interval, -3, PREV_NOT_SET);
      }
      if (max_p > end_2 && min_p <= end_2) {
        out.piece_list.emplace_back(SquareP, region_2 * LinearUP + LinearP,
                                    region_2 * region_2 * SquareU + region_2 * LinearU + Constant, end_2, max_p, -3, PREV_NOT_SET);
      }
      if (min_p > end_2) {
        out.piece_list.emplace_back(0, 0, INFINITY, min_p, max_p, PREV_NOT_SET, -3);
      }
      if (max_p < MACHINE_MAX) {
        out.piece_list.emplace_back(0, 0, INFINITY, max_p, MACHINE_MAX, PREV_NOT_SET, -3);
      }

  }


//crude checks

  double last_mean = MACHINE_MIN;
  for(SquareLossPieceList::iterator it = out.piece_list.begin(); it != out.piece_list.end(); it++) {

      if (it -> min_mean != last_mean) {
        printf("from original interval \n");
        print();
        printf("created this min interval\n");
        out.print();
        throw std::runtime_error("minimization error in constructing new piece. endpoint do not match");
      }

      if (it ->max_mean - it -> min_mean <= 0) {
        printf("from original interval \n");
        print();
        printf("created this min interval\n");
        out.print();
        throw std::runtime_error("minimization error in constructing new piece. segment length <= 0");
      }

      if (it -> Linear == -INFINITY || it -> Linear == INFINITY) {
        printf("from original interval \n");
        print();
        printf("created this min interval\n");
        out.print();
        throw std::runtime_error("minimization error in constructing new piece. linear coefs are Inf");
      }
      if (it -> Square == -INFINITY || it -> Square == INFINITY) {
        printf("from original interval \n");
        print();
        printf("created this min interval\n");
        out.print();
        throw std::runtime_error("minimization error in constructing new piece. square coefs are Inf");
      }

    last_mean = it -> max_mean;


  }

  return out;
}


void BiSquareLossPiece::print(){
  printf("%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",
         SquareU, LinearU, LinearUP, SquareP, LinearP, Constant,
         min_u, max_u, min_p, max_p);
}



void PiecewiseBiSquareLoss::print(){
  BiSquareLossPieceList::iterator it;
  printf("%10s %10s %15s %15s %15s %15s %15s %15s %15s %15s\n",
         "SquareU", "LinearU", "LinearUP", "SquareP", "LinearP", "Constant",
         "min_u", "max_u",
         "min_p", "max_p");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}

void PiecewiseBiSquareLoss::print_quad_p() {
  BiSquareLossPieceList::iterator it;
  printf("%10s %10s %15s %15s %15s %15s %s\n",
         "Square", "Linear", "Constant",
         "min_mean", "max_mean",
         "prev_mean", "data_i");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    printf("%.5f %.5f %.5f %.5f %.5f %.5f %d\n",
           it -> SquareP, it -> LinearP, it -> Constant,
           it -> min_p, it -> max_p,
           0.0,0);
  }
}

PiecewiseSquareLoss PiecewiseBiSquareLoss::get_univariate_p() {
  BiSquareLossPieceList::iterator it;
  PiecewiseSquareLoss out;
  out.piece_list.clear();

  for(it=piece_list.begin(); it != piece_list.end(); it++){
    out.piece_list.emplace_back(it -> SquareP, it -> LinearP, it -> Constant,
           it -> min_p, it -> max_p,
           0.0,0);
  }
  return out;
}

PiecewiseSquareLoss PiecewiseBiSquareLoss::min_over_u() {
  /*
   * min over each piece then 1d pruning
   *
   */
  int verbose = 0;
  BiSquareLossPieceList::iterator it = piece_list.begin();
  PiecewiseSquareLoss cost_min, cost_prev, cost_cur;
  int i = 0;

  while(it != piece_list.end()) {
    if (i == 0) {

      if (verbose) {
        printf("min of\n");
        printf("\t");
        it -> print();
      }


      cost_prev = it->min_over_u();

      if (verbose) {
        printf("min is .... \t");
        cost_prev.print();
      }


      cost_min = cost_prev;
    } else {

      if (verbose) {
        printf("min of\n");
        printf("\t");
        it -> print();
      }


      cost_cur = it->min_over_u();

      if (verbose) {
        printf("min is .... \t");
        cost_cur.print();

        printf("=cost current\n");
        cost_cur.print();
        printf("=cost previous\n");
        cost_prev.print();
      }



      cost_min.set_to_min_env_of(&cost_prev, &cost_cur, 0);
      int status = cost_min.check_min_of(&cost_prev, &cost_cur);
        if(status){
          cost_min.set_to_min_env_of(&cost_prev, &cost_cur, 1);
          printf("=cost current\n");
          cost_cur.print();
          printf("=cost previous\n");
          cost_prev.print();
          printf("=new cost model\n");
          cost_min.print();
          throw std::runtime_error("min_u error!");
        }
      cost_prev = cost_min;
    }
    i++;
    it++;
  }
  return cost_min;

}

void PiecewiseBiSquareLoss::set_to_pw_u(PiecewiseSquareLoss *f) {
  SquareLossPieceList::iterator it = f -> piece_list.begin();
  piece_list.clear();
  while(it != f -> piece_list.end()) {
    piece_list.emplace_back(it -> Square, it -> Linear, 0, 0, 0, it -> Constant, it -> min_mean, it -> max_mean, MACHINE_MIN, MACHINE_MAX);
    it++;
  }
}

void PiecewiseBiSquareLoss::set_to_pw_p(PiecewiseSquareLoss *f) {
  SquareLossPieceList::iterator it = f -> piece_list.begin();
  piece_list.clear();
  while(it != f -> piece_list.end()) {
    piece_list.emplace_back(0, 0, 0, it -> Square, it -> Linear, it -> Constant, MACHINE_MIN, MACHINE_MAX, it -> min_mean, it -> max_mean);
    it++;
  }
 }

void PiecewiseBiSquareLoss::add(double a, double b, double c, double d, double e, double f) {
  for(BiSquareLossPieceList::iterator it = piece_list.begin(); it != piece_list.end(); it++){
    it -> SquareU += a;
    it -> LinearU += b;
    it -> LinearUP += c;
    it -> SquareP += d;
    it -> LinearP += e;
    it -> Constant += f;
  }
}

void PiecewiseBiSquareLoss::set_to_addition_of(PiecewiseBiSquareLoss *f, PiecewiseBiSquareLoss *g, int verbose) {
  // fresh start
  piece_list.clear();

  BiSquareLossPieceList::iterator itf;
  BiSquareLossPieceList::iterator itg;

  // collect all unique interval endpoints in u and p
  std::set<double> u_pts, p_pts;
  for(itf = f->piece_list.begin(); itf != f->piece_list.end(); itf++) {
    u_pts.insert(itf -> min_u);
    u_pts.insert(itf -> max_u);
    p_pts.insert(itf -> min_p);
    p_pts.insert(itf -> max_p);
  }

  for(itg = g->piece_list.begin(); itg != g->piece_list.end(); itg++) {
    u_pts.insert(itg -> min_u);
    p_pts.insert(itg -> min_p);
    u_pts.insert(itg -> max_u);
    p_pts.insert(itg -> max_p);
  }

  std::set<double>::iterator it_u, it_p;
  int i = 0;
  int j = 0;
  double mu_min, mu_max, phi_min, phi_max;
  for (it_u=u_pts.begin(); it_u!=u_pts.end(); it_u++) {
    if (i == 0) {
      mu_min = *it_u;
      mu_max = *it_u;
    } else {
      mu_min = mu_max;
      mu_max = *it_u;
    }
    j = 0;
    for (it_p=p_pts.begin(); it_p!=p_pts.end(); it_p++) {
      if (j == 0 || i == 0) {
        phi_min = *it_p;
        phi_max = *it_p;
      } else {
        phi_min = phi_max;
        phi_max = *it_p;

        if (verbose) printf("searching matching functions for interval mu: [%f, %f] \t phi: [%f, %f] \n", mu_min, mu_max, phi_min, phi_max);
        itf = f->piece_list.begin();
        itg = g->piece_list.begin();
        while(!(itf-> min_u <= mu_min && itf->max_u >= mu_max && itf-> min_p <= phi_min && itf->max_p >= phi_max) && itf != f -> piece_list.end()) {
          if (verbose) {
            printf("\t search point failed at:\n");
            itf->print();
          }
          itf++;
        }

        while(!(itg-> min_u <= mu_min && itg->max_u >= mu_max && itg-> min_p <= phi_min && itg->max_p >= phi_max) && itg != g -> piece_list.end()) {
          if (verbose) {
            printf("\t search point failed at:\n");
            itg -> print();
          }
          itg++;
        }

        if ((itf-> min_u <= mu_min && itf->max_u >= mu_max && itf-> min_p <= phi_min && itf->max_p >= phi_max) &&
                (itg-> min_u <= mu_min && itg->max_u >= mu_max && itg-> min_p <= phi_min && itg->max_p >= phi_max)) { // success
          if (verbose) {
            printf("found matching functions for interval mu: [%f, %f] \t phi: [%f, %f] \n", mu_min, mu_max, phi_min, phi_max);
            itf->print();
            itg->print();
          }
          piece_list.emplace_back(itf -> SquareU + itg -> SquareU,
                                  itf -> LinearU + itg -> LinearU,
                                  itf -> LinearUP + itg -> LinearUP,
                                  itf -> SquareP + itg -> SquareP,
                                  itf -> LinearP + itg -> LinearP,
                                  itf -> Constant + itg -> Constant,
                                  mu_min, mu_max, phi_min, phi_max);
        }  else { // failed to find two functions over same domain
          if (verbose) {
            printf("NO matching functions for interval mu: [%f, %f] \t phi: [%f, %f] \n", mu_min, mu_max, phi_min, phi_max);
            f -> print();
            g -> print();
          }
          piece_list.emplace_back(0,
                                  0,
                                  0,
                                  0,
                                  0,
                                  INFINITY,
                                  mu_min, mu_max, phi_min, phi_max);
        }


      }
      j++;
    }
    i++;
  }


}



void PiecewiseBiSquareLosses::collect(PiecewiseBiSquareLoss * f) {
  piece_list.emplace_back(*f);
}

void PiecewiseBiSquareLosses::add(double a, double b, double c, double d, double e, double f) {
  for(PiecewiseBiSquareLossList::iterator it = piece_list.begin(); it != piece_list.end(); it++){
    it -> add(a, b, c, d, e, f);
  }
}

void PiecewiseBiSquareLosses::print() {
  printf("Printing whole collection\n");
  for(PiecewiseBiSquareLossList::iterator it = piece_list.begin(); it != piece_list.end(); it++){
    it -> print();
  }
}

PiecewiseBiSquareLoss PiecewiseBiSquareLosses::min_u() {
  int verbose = 0;
  PiecewiseBiSquareLossList::iterator it = piece_list.begin();
  PiecewiseSquareLoss cost_min, cost_prev, cost_cur;
  int i = 0;

  while(it != piece_list.end()) {
    if (i == 0) {
      cost_prev = it->min_over_u();
      cost_min = cost_prev;
    } else {
      cost_cur = it->min_over_u();
      cost_min.set_to_min_env_of(&cost_prev, &cost_cur, verbose);
      int status = cost_min.check_min_of(&cost_prev, &cost_cur);
      if(status){
          cost_min.set_to_min_env_of(&cost_prev, &cost_cur, 1);
          printf("=cost current\n");
          cost_cur.print();
          printf("=cost previous\n");
          cost_prev.print();
          printf("=new cost model\n");
          cost_min.print();
          throw std::runtime_error("PW-BiSquareLosses::min_u error!");
      }
      cost_prev = cost_min;
    }
    i++;
    it++;
  }
  PiecewiseBiSquareLoss out;
  out.set_to_pw_p(&cost_min);
  return out;
}


void PiecewiseBiSquareLosses::set_to_addition_of(PiecewiseBiSquareLosses *f, PiecewiseBiSquareLosses *g, int verbose) {
  // fresh start
  piece_list.clear();
  PiecewiseBiSquareLoss add;
  PiecewiseBiSquareLoss f_i, g_i;

  for(PiecewiseBiSquareLossList::iterator itf = f->piece_list.begin(); itf != f->piece_list.end(); itf++) {
    f_i = *itf;
    for(PiecewiseBiSquareLossList::iterator itg = g->piece_list.begin(); itg != g->piece_list.end(); itg++) {

      g_i = *itg;

      if (verbose) {
        printf("--------\n");
        printf("current f_i\n");
        f_i.print();
        printf("current g_i\n");
        g_i.print();
      }

      add.set_to_addition_of(&f_i, &g_i, verbose);

      if (verbose) {
        printf("addition\n");
        add.print();
      }

      piece_list.emplace_back(add);
    }
  }
}