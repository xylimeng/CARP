// Copyright (C) 2014 National ICT Australia (NICTA)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
//
// Written by Conrad Sanderson - http://conradsanderson.id.au
// Written by George Yammine


// Demonstration of how to connect Armadillo with Matlab mex functions.
// Version 0.2


#include "armaMex.hpp"
#include "helper.hpp"
#include "tree_class.hpp"

using namespace std;
using namespace arma;


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check the number of input arguments.
    if (nrhs != 3)
        mexErrMsgTxt("Incorrect number of input arguments.");
    
    vec obs = armaGetPr(prhs[0]);
    vec dimension_vec = armaGetPr(prhs[1]);
    uvec dimension = conv_to<uvec>::from( dimension_vec );
    
    vec hyper = armaGetPr(prhs[2]);  
            
    // Perform calculations
    class_tree tree_test(obs, dimension);
    mat par = tree_test.hyper2par(hyper); 
    vec log_lambda_d = tree_test.get_lambda_map(par); 
    vec log_M_d = tree_test.get_M_d_map(par); 
    
    uvec node_idx_d = tree_test.node_idx_d; 
    uvec left_idx_d = tree_test.left_idx_d; 
    uvec right_idx_d = tree_test.right_idx_d; 
    umat C(node_idx_d.n_elem, 3); 
    C.col(0) = node_idx_d; 
    C.col(1) = left_idx_d; 
    C.col(2) = right_idx_d; 
    
    uvec compact_idx_d = tree_test.compact_idx_map; 
    uvec flag_first = tree_test.flag_first(0); 
    uvec flag_not_first = tree_test.flag_not_first(0); 
    umat pointer_flag(2, (tree_test.flag_first).n_elem); // cumulative lengths; 
    pointer_flag(0, 0) = flag_first.n_elem; 
    pointer_flag(1, 0) = flag_not_first.n_elem;
    
    for (int i = 1; i < (tree_test.flag_first).n_elem; i++){
        flag_first.insert_rows(pointer_flag(0, i - 1), tree_test.flag_first(i)); 
        flag_not_first.insert_rows(pointer_flag(1, i - 1), tree_test.flag_not_first(i)); 
        pointer_flag(0, i) = flag_first.n_elem; 
        pointer_flag(1, i) = flag_not_first.n_elem;
    }
    
    umat rank_left_child = tree_test.rank_left_child_2; 
    umat rank_right_child = tree_test.rank_right_child_2; 
    umat dividible = tree_test.dividible_2; 
    
    vec y = tree_test.y; 
    vec ss = tree_test.ss; 
    
  
    
            
    // Create the output argument plhs[0] to return cube C
    plhs[0] = armaCreateMxMatrix(log_lambda_d.n_rows, log_lambda_d.n_cols);
    plhs[1] = armaCreateMxMatrix(C.n_rows, C.n_cols, mxUINT64_CLASS); 
    plhs[2] = armaCreateMxMatrix(compact_idx_d.n_elem, 1, mxUINT64_CLASS); 
    plhs[3] = armaCreateMxMatrix(flag_first.n_elem, 1, mxUINT64_CLASS); 
    plhs[4] = armaCreateMxMatrix(flag_not_first.n_elem, 1, mxUINT64_CLASS); 
    plhs[5] = armaCreateMxMatrix(pointer_flag.n_rows, pointer_flag.n_cols, mxUINT64_CLASS); 
    plhs[6] = armaCreateMxMatrix(rank_left_child.n_rows, rank_left_child.n_cols, mxUINT64_CLASS); 
    plhs[7] = armaCreateMxMatrix(rank_right_child.n_rows, rank_right_child.n_cols, mxUINT64_CLASS); 
    plhs[8] = armaCreateMxMatrix(dividible.n_rows, dividible.n_cols, mxUINT64_CLASS); 
    plhs[9] = armaCreateMxMatrix(par.n_rows, par.n_cols);
    plhs[10] = armaCreateMxMatrix(y.n_rows, y.n_cols);
    plhs[11] = armaCreateMxMatrix(ss.n_rows, ss.n_cols);
    plhs[12] = armaCreateMxMatrix(log_M_d.n_rows, log_M_d.n_cols);
    
    // Return the cube C as plhs[0] in Matlab/Octave
    armaSetPr(plhs[0], log_lambda_d);
    armaSetData(plhs[1], C); // non-double data uses "armaSetData" 
    armaSetData(plhs[2], compact_idx_d); 
    armaSetData(plhs[3], flag_first); 
    armaSetData(plhs[4], flag_not_first); 
    armaSetData(plhs[5], pointer_flag); 
    armaSetData(plhs[6], rank_left_child); 
    armaSetData(plhs[7], rank_right_child); 
    armaSetData(plhs[8], dividible); 
    armaSetPr(plhs[9], par);
    armaSetData(plhs[10], y); 
    armaSetData(plhs[11], ss); 
    armaSetData(plhs[12], log_M_d);
    
    return;
}
