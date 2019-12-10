

#include "armaMex.hpp"
#include "helper.hpp"
#include "tree_class.hpp"

using namespace std;
using namespace arma;


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check the number of input arguments.
    if (nrhs != 4)
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
    
    mat par = armaGetPr(prhs[2]);
    int step = armaGetScalar<int>(prhs[3]);
    
    // Perform calculations
    class_tree tree_test(obs, dimension);
    mat BMA = tree_test.fit_tree_cs(par, step);
    
    // Create the output argument plhs[0] to return cube C
    plhs[0] = armaCreateMxMatrix(BMA.n_rows, BMA.n_cols);
    
    // Return the cube C as plhs[0] in Matlab/Octave
    armaSetPr(plhs[0], BMA);
    
    return;
}
