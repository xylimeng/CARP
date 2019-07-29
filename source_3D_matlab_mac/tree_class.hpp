//
//  tree_class.hpp
//  Surfing
//
//  Created by Meng Li on 7/15/16.
//  Copyright © 2016 Meng Li. All rights reserved.
//

#ifndef tree_class_hpp
#define tree_class_hpp

// #include <mex.h>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include </usr/local/include/armadillo>
#include <random>

using namespace std;
using namespace arma;

// define the tree class
class class_tree
{
public:
    vec obs;
    uvec dimension;
    
    class_tree(const vec& obs_input, uvec& dimension_input); // constructor: generate dictionary
    
    vec get_phi_map(const mat& par, const vec& w_d, const vec& ss_adjust); // obtain log_phi
    
    vec get_M_d_map(const mat& par); // obtain log(M_d)
    
    // vec get_lambda_map(const mat& par); // obtain log(lambda_d)
    void get_post_map(const mat& par); // obtain log(lambda_d) and post_eta (log scale); both posterior maps
    
    double marginal_likelihood(const mat& par); // marginal likelihood
    
    mat hyper2par(const vec &x); // re-parameterization for optim
    
    vec fit_tree(const mat& par, const vec& shift_vec); // estimtation without cycle spinning
    vec fit_tree(const mat& par); // estimtation without cycle spinning; default: no shift 
    vec fit_tree_cs(const mat& par, int step); // estimation via cycle spinning;
    
//private:
    uword N, total_level, m, pair_N, N_prime;
    uvec node_idx_d, left_idx_d, right_idx_d, compact_idx_map, num_child, level, pair_level, pair_direction;
    field<uvec> flag_first, flag_not_first, F;
    vec n_A, num_child_double, ss_adjust_orig, w_d_orig;
    vec y, ss; 
    
    // used for matlab trials - didn't export these var's in the past 
    umat rank_left_child_2, rank_right_child_2, dividible_2; //'_2': questionable - remove finally
    umat rank_left_child, rank_right_child, dividible; 
    double rescale_sum;
    
    // draw trees
    // dictionary for family labels (or 'rank'): 'rank' starts from 1
    umat family_rank;
    uvec position; // position in the raw location space; useful for the last scale: starts from 0
    vec post_lambda_d;
    vec post_eta; 
    mat lambda_mat;
    void get_lambda_mat(); // obtain matrix form of lambda map (a lot of zeros; not compact)
//    umat draw_post_position(const uword& n_smp); // draw samples from the posterior distribution of tree; return position
    umat draw_post_position(const uword& n_smp); // draw samples from the posterior distribution of tree; return (direction, pruning indicator, position) by (number of samples)
    
    umat MAP_tree; mat MAP_fit;
    void fit_MAP_tree(const mat& par, const vec& shift_vec); // obtain the MAP tree and estimate
    void fit_MAP_tree(const mat& par); // obtain the MAP tree and estimate; no shift as default
};

#endif /* tree_class_hpp */
