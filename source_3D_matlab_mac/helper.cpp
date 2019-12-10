

#include "helper.hpp"

// const double log_inv_sqrt_2pi = std::log( 1 / sqrt(2*M_PI) );

// [[Rcpp::export]]
double normal_logpdf(double x, double s)
{
    double a = x / s;
    return log_inv_sqrt_2pi - log(s) - 0.5 * a * a;
}

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

// [[Rcpp::export]]
vec log_exp_x_plus_exp_y_vec(const arma::vec& c, const arma::vec& d){
    vec x(c.n_elem);
    for (int i = 0; i < c.n_elem; i++)
    {
        x(i) = log_exp_x_plus_exp_y(c(i), d(i));
    }
    return x;
    
}

uvec merge_left_right (const uvec& v1, const uvec& v2)
{
    uvec v3(v1.n_elem + v2.n_elem);
    uvec v4 = regspace<uvec>(0,  2,  v3.n_elem - 1);
    v3(v4) = v1; v3(v4 + 1) = v2;
    return v3;
};

vec merge_left_right (const vec& v1, const vec& v2)
{
    vec v3(v1.n_elem + v2.n_elem);
    uvec v4 = regspace<uvec>(0,  2,  v3.n_elem - 1);
    v3(v4) = v1; v3(v4 + 1) = v2;
    return v3;
};

vec mu_one(vec x, double scale)
{
    // Gaussian: scale = tau.j
    return x / (1 + 1/scale);
    
}

vec mu_one(vec x, vec scale)
{
    // Gaussian: scale = tau.j
    return x / (1 + 1/scale);
    
}

// rank: id of each node; won't be useful otherwise
// idx: idx to find nodes. Depend on platform. For c++, idx starts from 0; for R or matlab, idx starts from 1.
//      for c++: idx = rank - 1; for R/Matlab, idx = rank.

// [[Rcpp::export]]
umat rank_vec2tube(const uvec& rank, const uvec& num_node)
{
    uword m = num_node.n_elem;
    umat tube_rank(rank.n_elem, m);
    
    // easier to start everything from 0 - local change
    // initialize tube_rank: map from t to s = (s_1, ..., s_m): t = rank
    uvec current_rank = rank - 1;
    
    for (int d = 0; d < m; d++) // use linear indexing to transfer 'rank' to 'tube_rank'
    {
        uvec residule = floor(current_rank / num_node(d));
        tube_rank.col(d) = current_rank - residule * num_node(d); // obtain s_d
        
        current_rank = residule;
    }
    
    tube_rank = tube_rank + 1;
    return(tube_rank);
}

// [[Rcpp::export]]
// recover the rank again
uvec rank_tube2vec(const umat& tube_rank, const uvec& num_node)
{
    uword m = num_node.n_elem;
    uvec recover_rank = tube_rank.col(0) - 1;
    uvec temp = cumprod(num_node);
    for (int d = 1; d < m; d++)
    {
        recover_rank = recover_rank + (tube_rank.col(d) - 1) * temp(d - 1);
    }
    recover_rank = recover_rank + 1;
    return(recover_rank);
}

// [[Rcpp::export]]
double log_ratio_function(double x, double scale){
    double ret = log(1/sqrt(1 + scale)) + (0.5 * (x * x) / (1 + 1/scale));
    return ret;
}

double MSE(const mat &x, const mat &y)
{
    return accu((x - y) % (x - y))/x.n_elem * 100;
}

mat circshift(const mat&xx, const vec&shift)
{
    mat ret = xx;
    for (int d=0; d<shift.n_elem; d++)
    {
        ret = arma::shift(ret, shift(d), d);
    };
    return ret;
}



umat ind2sub(const uvec& ind, const uvec& dimension)
{
    uword m = dimension.n_elem;
    umat sub(ind.n_elem, m);
    
    // easier to start everything from 0 - local change
    // initialize sub: map from t to s = (s_1, ..., s_m): t = ind
    
    uvec ind_copy = ind;
    
    for (int d = 0; d < m; d++) // use linear indexing to transfer 'ind' to 'sub'
    {
        for(int i = 0; i < ind.n_elem; ++i){
            sub(i, d) = ind_copy(i) % dimension(d); // obtain s_d
 //           ind(i) = ind(i) / dimension(d);
        }
        
        ind_copy /= dimension(d);
    }
    
    return(sub);
}

uvec ind2sub(const uword& ind, const uvec& dimension)
{
    uword m = dimension.n_elem;
    uvec sub(m);
    
    // easier to start everything from 0 - local change
    // initialize sub: map from t to s = (s_1, ..., s_m): t = ind
    
    uword ind_copy = ind;
    
    for (int d = 0; d < m; d++) // use linear indexing to transfer 'ind' to 'sub'
    {

            sub(d) = ind_copy % dimension(d); // obtain s_d
            //           ind(i) = ind(i) / dimension(d);
        
        ind_copy /= dimension(d);
    }
    
    return(sub);
}

uvec sub2ind(const umat& sub, const uvec& dimension)
{
    uword m = dimension.n_elem;
    uvec ind = sub.col(0);
    uvec temp = cumprod(dimension);
    for (int d = 1; d < m; d++)
    {
        ind += sub.col(d) * temp(d - 1);
    }
    return(ind);
}

vec circshift_vectorize(const vec& X, const uvec& dimension, const vec& shift_vec)
{
    uword num_elem = X.n_elem;
    uword m = dimension.n_elem;
    
    uvec ind = linspace<uvec>(0, num_elem - 1, num_elem);
    // step 1: return the subscript index ('ind2sub')
    umat sub = ind2sub(ind, dimension);
    // step 2: performe circulant shift
    for (int d = 0; d < dimension.n_elem; ++d)
    {
        sub.col(d) +=  (dimension(d) - shift_vec(d)); // or sub.col(d) -= shift_vec(d); but maybe contain bugs by using minus so add dimension(d);
        for (int i = 0; i < num_elem; ++i){
            sub(i, d) %= dimension(d);
        }
    }
    // step 3: return the linear index ('sub2ind')
    uvec shift_ind = sub2ind(sub, dimension);
    vec X_shift = X(shift_ind);
    
    return(X_shift);
}
