


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
    uvec rank = conv_to<uvec>::from( armaGetPr(prhs[0]) );  
    uvec num_node = conv_to<uvec>::from( armaGetPr(prhs[1]) ); 
           
    
    umat z = rank_vec2tube(rank, num_node);
    
            
    // Create the output argument plhs[0] to return cube C
    plhs[0] = armaCreateMxMatrix(z.n_rows, z.n_cols, mxUINT64_CLASS);

    
    // Return the cube C as plhs[0] in Matlab/Octave
    armaSetData(plhs[0], z);
 
    return;
}
