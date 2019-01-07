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
    if (nrhs != 2)
        mexErrMsgTxt("Incorrect number of input arguments.");
    
    // Check type of input.
    if ( (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS))
        mexErrMsgTxt("Input must me of type double.");
    
    // Check if input is real.
    if ( (mxIsComplex(prhs[0])))
        mexErrMsgTxt("Input must be real.");
    
    // Create matrices X and Y from the first and second argument. 
    vec x = armaGetPr(prhs[0]);  
    vec y = armaGetPr(prhs[1]); 
    
    vec z = log_exp_x_plus_exp_y_vec(x, y); 
           
    
  
    
            
    // Create the output argument plhs[0] to return cube C
    plhs[0] = armaCreateMxMatrix(z.n_rows, z.n_cols);

    
    // Return the cube C as plhs[0] in Matlab/Octave
    armaSetPr(plhs[0], z);
 
    return;
}
