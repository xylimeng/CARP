//
//  helper.hpp
//  Surfing
//
//  Created by Meng Li on 7/15/16.
//  Copyright © 2016 Meng Li. All rights reserved.
//

#ifndef helper_hpp
#define helper_hpp

// #include <mex.h>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <math.h>
#include </dscrhome/ml371/privatemodules/usr/include/armadillo>

using namespace std;
using namespace arma;

const double log_inv_sqrt_2pi = std::log( 1 / sqrt(2*M_PI) );

double normal_logpdf(double x, double s);

double log_exp_x_plus_exp_y(double x, double y);

vec log_exp_x_plus_exp_y_vec(const arma::vec& c, const arma::vec& d);

uvec merge_left_right (const uvec& v1, const uvec& v2);

vec merge_left_right (const vec& v1, const vec& v2);

vec mu_one(vec x, double scale);

vec mu_one(vec x, vec scale);

umat rank_vec2tube(const uvec& rank, const uvec& num_node);
// recover the rank again
uvec rank_tube2vec(const umat& tube_rank, const uvec& num_node);

double log_ratio_function(double x, double scale);

double MSE(const mat &x, const mat &y);

mat circshift(const mat&xx, const vec&shift);

umat ind2sub(const uvec& ind, const uvec& dimension);
uvec sub2ind(const umat& sub, const uvec& dimension);
uvec ind2sub(const uword& ind, const uvec& dimension);
vec circshift_vectorize(const vec& X, const uvec& dimension, const vec& shift_vec);  // returns circshifted X: both input and output are vectorized using linear indexing.



#endif /* helper_hpp */
