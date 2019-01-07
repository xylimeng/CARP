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
    
    // Check type of input.
    if ( (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS))
        mexErrMsgTxt("Input must me of type double.");
    
    // Check if input is real.
    if ( (mxIsComplex(prhs[0])))
        mexErrMsgTxt("Input must be real.");
    
    // Create matrices X and Y from the first and second argument.
    vec obs = armaGetPr(prhs[0]);
    vec dimension_vec = armaGetPr(prhs[1]);
    uvec dimension = conv_to<uvec>::from( dimension_vec );
    mat hyper_grid = armaGetPr(prhs[2]);
    
    
    // Perform calculations
    class_tree tree_test(obs, dimension);
    
    int n_grid = hyper_grid.n_cols;
    vec ml_mat(n_grid);

    for (int i = 0; i < n_grid; i++)
    {
        vec hyper = hyper_grid.col(i);
        mat par = tree_test.hyper2par(hyper);
        double ml = tree_test.marginal_likelihood(par);
        ml_mat(i) = ml;
    }
    
    
    // Create the output argument plhs[0] to return cube C
    plhs[0] = armaCreateMxMatrix(ml_mat.n_rows, ml_mat.n_cols);
    
    // Return the cube C as plhs[0] in Matlab/Octave
    armaSetPr(plhs[0], ml_mat);
    
    return;
}
