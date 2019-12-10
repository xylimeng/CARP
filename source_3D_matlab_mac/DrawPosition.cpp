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
    
    arma_rng::set_seed(seed); //set_seed 2019;
    // arma_rng::set_seed_random() // set random seed
    
    umat post_position = tree_test.draw_post_position(n_smp);
    
    // Create the output argument plhs[0] to return cube C
    plhs[0] = armaCreateMxMatrix(post_position.n_rows, post_position.n_cols, mxUINT64_CLASS);
    
    // Return the cube C as plhs[0] in Matlab/Octave
    armaSetData(plhs[0], post_position);
    
    return;
}
