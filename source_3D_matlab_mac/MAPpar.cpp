


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
    
    mat par = armaGetPr(prhs[2]);
            
    // Perform calculations
    class_tree tree_test(obs, dimension);
    
    tree_test.get_post_map(par);
    tree_test.get_lambda_mat();
    
    tree_test.fit_MAP_tree(par);
    
    umat MAP_tree = tree_test.MAP_tree;
    vec MAP_fit = tree_test.MAP_fit;
    
    // Create the output argument plhs[0] to return cube C
    plhs[0] = armaCreateMxMatrix(MAP_tree.n_rows, MAP_tree.n_cols, mxUINT64_CLASS);
    plhs[1] = armaCreateMxMatrix(MAP_fit.n_rows, MAP_fit.n_cols);
    
    // Return the cube C as plhs[0] in Matlab/Octave
    armaSetData(plhs[0], MAP_tree);
    armaSetData(plhs[1], MAP_fit);
    
    return;
}
