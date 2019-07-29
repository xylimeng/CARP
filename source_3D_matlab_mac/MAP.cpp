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
    if (nrhs != 5)
        mexErrMsgTxt("Incorrect number of input arguments.");
    
    vec obs = armaGetPr(prhs[0]);
    vec dimension_vec = armaGetPr(prhs[1]);
    uvec dimension = conv_to<uvec>::from( dimension_vec );
    
    vec hyper = armaGetPr(prhs[2]);
    int n_smp = armaGetScalar<int>(prhs[3]);
    int seed = armaGetScalar<int>(prhs[4]);
            
    // Perform calculations
    class_tree tree_test(obs, dimension);
    mat par = tree_test.hyper2par(hyper);
    
    tree_test.get_post_map(par);
    tree_test.get_lambda_mat();
    
    tree_test.fit_MAP_tree(par);
    
    umat MAP_tree = tree_test.MAP_tree;
    vec MAP_fit = tree_test.MAP_fit;
    
    // Create the output argument plhs[0] to return cube C
    plhs[0] = armaCreateMxMatrix(MAP_tree.n_rows, MAP_tree.n_cols, mxUINT64_CLASS);
    plhs[1] = armaCreateMxMatrix(MAP_tree.n_rows, MAP_tree.n_cols);
    
    // Return the cube C as plhs[0] in Matlab/Octave
    armaSetData(plhs[0], MAP_tree);
    armaSetData(plhs[1], MAP_fit);
    
    return;
}
