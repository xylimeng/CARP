#include <cmath>
#include <math.h>
#include "matrix.h"  // mwSize, mwIndex
#include <mex.h>
#include <cstring>

using namespace std;

double log_exp_x_plus_exp_y(double x, double y)
{
    double result;
    if( ( std::isinf( fabs(x) ) == true ) && ( std::isinf( fabs(y) ) == false )  )
        result = y;
    else if ( ( std::isinf( fabs(x) ) == false ) && ( std::isinf( fabs(y) ) == true )  )
        result = x;
    else if ( ( std::isinf( fabs(x) ) == true ) && ( std::isinf( fabs(y) ) == true )  )
        result = x;
    else if ( x - y >= 100 ) result = x;
    else if ( x - y <= -100 ) result = y;
    else {
        if (x > y) {
            result = y + log( 1 + exp(x-y) );
        }
        else result = x + log( 1 + exp(y-x) );
    }
    return result;
}

void mexFunction(int nlhs, mxArray *plhs[],  // plhs: array of mxArray pointers to output
        int nrhs, const mxArray *prhs[])    // prhs: array of mxArray pointers to input
{
    const mwSize *dims;
    dims = mxGetDimensions(prhs[0]);
    size_t elements=mxGetNumberOfElements(prhs[0]);
    
    // Macros for the ouput and input arguments
    
    double  *pointer;          /* pointer to real data in new array */
    double *pointer_x;       // pointers to input matrix
    double *pointer_y;
    
    /* Create an m-by-n mxArray; you will copy existing data into it */
    plhs[0] = mxCreateNumericMatrix(dims[0], dims[1], mxDOUBLE_CLASS, mxREAL);
    pointer = mxGetPr(plhs[0]);
    
    pointer_x = (double *)mxGetData(prhs[0]);
    pointer_y = (double *)mxGetData(prhs[1]);

    for(int j=0; j<elements; j++){
        pointer[j] = log_exp_x_plus_exp_y(pointer_x[j], pointer_y[j]);
    }
    
    
    return;
}