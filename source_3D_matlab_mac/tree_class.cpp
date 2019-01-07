//
//  tree_class.cpp
//  Surfing
//
//  Created by Meng Li on 7/15/16.
//  Copyright © 2016 Meng Li. All rights reserved.
//

#include "tree_class.hpp"
#include "helper.hpp"


class_tree::class_tree (const vec& obs_input, uvec& dimension_input)
{
    obs = obs_input;
    dimension = dimension_input; 
    Col<uword> J_vec = log2(dimension);
    uvec num_node = 2 * dimension - 1; //total number of nodes
    
    N = prod(num_node);
    m = J_vec.n_elem;
     
    uvec rank = linspace<uvec>(1, N, N);
    umat tube_rank = rank_vec2tube(rank, num_node); //intermediate step, will not output
    umat tube_level = floor(log2(tube_rank)); //intermediate step, will not output
    Mat<uword> rank_left_child(N, m), rank_right_child(N, m), rank_parent(N, m), dividible(N, m), root(N, m);
    umat left_child = tube_rank - 2 * floor(tube_rank / 2); // (tube_rank % 2 == 0)
    
    for (int d=0; d < m; d++)
    {
        tube_rank.col(d) *= 2;
        rank_left_child.col(d) = rank_tube2vec(tube_rank, num_node);
        tube_rank.col(d) += 1;
        rank_right_child.col(d) = rank_tube2vec(tube_rank, num_node);
        tube_rank.col(d) -= 1;
        tube_rank.col(d) /= 2;
        
        umat temp_temp = tube_rank;
        temp_temp.col(d) /= 2;
        rank_parent.col(d) = rank_tube2vec(floor(temp_temp), num_node);
        dividible.col(d) = (tube_level.col(d) < J_vec(d));
        
        root.col(d) = (tube_rank.col(d) == 0);
    }
    
    rank_left_child_2 = rank_left_child; 
    rank_right_child_2 = rank_right_child; 
    dividible_2 = dividible; 
    
    level = sum(tube_level, 1); // calculate grand level
    
    num_child = sum(dividible, 1);
    
    //    uvec level = sum(tube_level, 1); // calculate grand level
    //    vec num_child = conv_to< vec >::from(sum(dividible, 1));
    
    total_level = sum(J_vec);
    
    // each level: return a field
    F.set_size(total_level + 1);
    // field<uvec> F(total_level + 1); // store level 0 to 'total_level'
    for (int l = 0; l <= total_level; l++){
        F(l) = (find(level == l));
    };
    
    // information for 'families': nodes that have (parents, children)
    // initialize it at d = 0
    // first generate N by d tube, then extract the 'temp' elements
    // for high dimension: use loops over dimension and 'push.back' instead of expanding to tubes
    uvec temp = find(dividible > 0);
    umat temp_node_idx_d = repmat(rank - 1, 1, m);
    node_idx_d = temp_node_idx_d(temp);
    
    left_idx_d = rank_left_child(temp) - 1;
    right_idx_d = rank_right_child(temp) - 1;
    pair_N = temp.n_elem;
    umat temp_direction = repmat(linspace<uvec>(0, m - 1, m), 1, N);
    umat temp_directioin_2 = temp_direction.t();
    pair_direction = temp_directioin_2(temp);
    umat rep_level = repmat(level, 1, m);
    pair_level = rep_level(temp);
    //    mat rep_num_child = repmat(num_child, 1, m);
    //    vec pair_num_child = rep_num_child(temp);
    umat cum_sum_dividible = cumsum(dividible, 1); //cumsum of each row
    uvec pair_cumsum_direction = cum_sum_dividible(temp); // temporary need - clean finally
    uvec pair_dividible = dividible(temp);
    
    // calculate F(l), flag_first, flag_not_first only for non-atom nodes
    // use atom-nodes for initilization
    // field<uvec> flag_first(total_level + 1, m + 1); // dth direction & first cut
    // field<uvec> flag_not_first(total_level + 1, m);  // dth direction & not first cut
    flag_first.set_size((total_level + 1) * (m + 1));
    flag_not_first.set_size((total_level + 1) * (m + 1));
    // R <-> cpp: transfer: as<field>: there is a bug
    // such that flag_first becomes a field with 1 column
    // we vectorize the index
    // (x, y) -> x + y * (total_level + 1)
    for (int l = 0; l <= total_level; l++){
        for (int d = 0; d < m; d++)
        {
            flag_first(l + d * (total_level + 1)) = find((pair_level == l) % (pair_direction == d) % (pair_cumsum_direction == 1));
            flag_not_first(l + d * (total_level + 1)) = find((pair_level == l) % (pair_direction == d) % (pair_cumsum_direction > 1));
            
        };
        
        //        uvec flag_first_temp = flag_first(l);
        //        for (int d = 1; d < m; d++)
        //        {
        //            flag_first_temp.insert_rows(flag_first_temp.n_elem, flag_first(l + d * (total_level + 1)));
        //        };
        flag_first(l + m * (total_level + 1)) = find((pair_level == l)  % (pair_cumsum_direction == 1));
        flag_not_first(l + m * (total_level + 1)) = find(pair_level == l);
        
    };
    
    // 'rank' of node is not sorted by 'level'; we need to create a map to find the index of parent nodes
    uvec atom_indicator = (num_child > 0);
    compact_idx_map = cumsum(atom_indicator) - 1;
    compact_idx_map(find(num_child == 0)).fill(-1); // '-1' is not supposed to be used; put here to detect bugs
    
    // all nodes: N; nodes that are d-type parent(family) - pair_N; nodes that are parent - N_prime;
    // pair_N is more than N_prime: in pair_N, nodes replicated if they are different types of parent
    
    N_prime = N - prod(dimension); // number of nodes except the last level
    vec v1 = arma::conv_to< vec >::from(total_level - level);
    n_A = exp2(v1);
    // v = sqrt(n_A);
    
    num_child_double = arma::conv_to< vec >::from(num_child);
    
    // rescale sum
    rescale_sum = accu(obs) / sqrt(obs.n_elem);
    // calculate 'obs' (y), 'ss' (ss), 'w_d' (w) for each node - bottom up
    y.set_size(N); ss.set_size(N); 
    // vec y(N, fill::zeros), ss(N, fill::zeros);
    y(F(total_level)) = obs;
    ss(F(total_level)) = obs % obs;
    
    for (int l = total_level - 1; l >= 0;l --){ // current level
        uvec index = flag_first(l + m * (total_level + 1)); // select the current node has dth cut & first cut
        uvec left_index = left_idx_d(index);
        uvec right_index = right_idx_d(index);
        uvec idx_node = node_idx_d(index);
        y(idx_node) = y(left_index) + y(right_index);
        ss(idx_node) = ss(left_index) + ss(right_index);
    }
    
    ss_adjust_orig = ss - ((y % y) / n_A);
    w_d_orig = y(left_idx_d) - y(right_idx_d);
    w_d_orig.each_col() /= sqrt(n_A(node_idx_d));
    
}


vec class_tree::get_phi_map(const mat& par, const vec& w_d, const vec& ss_adjust)
{
    // parameters: 3 rows by (total_level + 1): last element is sigma
    vec rho_vec = par.col(0);
    vec tau_vec = par.col(1);
    vec eta_vec = par.col(2);
    double sigma = eta_vec(total_level);
    eta_vec.shed_row(total_level);
    
    vec M_d(pair_N);
    vec log_rho = log(rho_vec);
    vec log_rho_c = log(1 - rho_vec);
    
    for (int i = 0; i < pair_N; i++)
    {
        M_d(i) = log_exp_x_plus_exp_y(log_rho(pair_level(i)) + log_ratio_function(w_d(i) / sigma, tau_vec(pair_level(i))),
                                      log_rho_c(pair_level(i)) ) + normal_logpdf(w_d(i), sigma);
        
    };
    
    vec log_Phi(N_prime); // from node 0 to node (N_prime - 1);
    vec log_p0 = (n_A - 1) * (log_inv_sqrt_2pi - log(sigma)) - (0.5 * (1/sigma) * (1/sigma)) * ss_adjust;
    
    // eta_vec(total_level - 1) not used at all
    
    for (int l = total_level - 1; l >= 0; l --) // current level
    {
        double eta = eta_vec(l);
        
        if (l == (total_level - 1))
            // initialize log_Phi: level (J - 1)
            // features: children levels log_Phi = 0;
            //  nodes (current level) only has 1 cut direction -> lambda_d = 1, num_child = 1
        {
            for (int d = 0; d < m; d++)
            {
                uvec idx_d = flag_first(l + d * (total_level + 1)); // if the current node has dth cut & first cut
                uvec idx_node = compact_idx_map(node_idx_d(idx_d));
                log_Phi(idx_node) = M_d(idx_d);
            };
        }
        else
        {
            // update log_Phi
            for (int d = 0; d < m; d++)
            {
                
                uvec index = flag_first(l + d * (total_level + 1)); // if the current node has dth cut & first cut
                uvec left_index = compact_idx_map(left_idx_d(index));
                uvec right_index = compact_idx_map(right_idx_d(index));
                uvec idx_node = compact_idx_map(node_idx_d(index));
                log_Phi(idx_node) = log_Phi(left_index) + log_Phi(right_index) + M_d(index);
                
                index = flag_not_first(l + d * (total_level + 1)); // if the current node has more than dth cut
                left_index = compact_idx_map(left_idx_d(index));
                right_index = compact_idx_map(right_idx_d(index));
                idx_node = compact_idx_map(node_idx_d(index));
                log_Phi(idx_node) = log_exp_x_plus_exp_y_vec(log_Phi(idx_node),
                                                             log_Phi(left_index) + log_Phi(right_index) + M_d(index));
            };
        }
        
        uvec idx_node = F(l);
        uvec compact_idx_node = compact_idx_map(idx_node);
        log_Phi(compact_idx_node) -= log(num_child_double(idx_node));
        
        log_Phi(compact_idx_node) = log_exp_x_plus_exp_y_vec(log(eta) + log_p0(idx_node), log(1 - eta) + log_Phi(compact_idx_node));
        
    };
    
    return log_Phi;
};

// [[Rcpp::export]]
double class_tree::marginal_likelihood(const mat& par)
{
    // parameters: 3 rows by (total_level + 1): last element is sigma
    vec rho_vec = par.col(0);
    vec tau_vec = par.col(1);
    vec eta_vec = par.col(2);
    double sigma = eta_vec(total_level);
    eta_vec.shed_row(total_level);
    
    vec M_d(pair_N);
    vec log_rho = log(rho_vec);
    vec log_rho_c = log(1 - rho_vec);
    
    vec w_d(N, fill::zeros), ss_adjust(N, fill::zeros);
    w_d = w_d_orig;
    ss_adjust = ss_adjust_orig;
    
    for (int i = 0; i < pair_N; i++)
    {
        M_d(i) = log_exp_x_plus_exp_y(log_rho(pair_level(i)) + log_ratio_function(w_d(i) / sigma, tau_vec(pair_level(i))),
                                      log_rho_c(pair_level(i)) ) + normal_logpdf(w_d(i), sigma);
        
    };
    
    vec log_Phi(N_prime); // from node 0 to node (N_prime - 1);
    vec log_p0 = (n_A - 1) * (log_inv_sqrt_2pi - log(sigma)) - (0.5 * (1/sigma) * (1/sigma)) * ss_adjust;
    
    // eta_vec(total_level - 1) not used at all
    
    for (int l = total_level - 1; l >= 0; l --) // current level
    {
        double eta = eta_vec(l);
        
        if (l == (total_level - 1))
            // initialize log_Phi: level (J - 1)
            // features: children levels log_Phi = 0;
            //  nodes (current level) only has 1 cut direction -> lambda_d = 1, num_child = 1
        {
            for (int d = 0; d < m; d++)
            {
                uvec idx_d = flag_first(l + d * (total_level + 1)); // if the current node has dth cut & first cut
                uvec idx_node = compact_idx_map(node_idx_d(idx_d));
                log_Phi(idx_node) = M_d(idx_d);
            };
        }
        else
        {
            // update log_Phi
            for (int d = 0; d < m; d++)
            {
                
                uvec index = flag_first(l + d * (total_level + 1)); // if the current node has dth cut & first cut
                uvec left_index = compact_idx_map(left_idx_d(index));
                uvec right_index = compact_idx_map(right_idx_d(index));
                uvec idx_node = compact_idx_map(node_idx_d(index));
                log_Phi(idx_node) = log_Phi(left_index) + log_Phi(right_index) + M_d(index);
                
                index = flag_not_first(l + d * (total_level + 1)); // if the current node has more than dth cut
                left_index = compact_idx_map(left_idx_d(index));
                right_index = compact_idx_map(right_idx_d(index));
                idx_node = compact_idx_map(node_idx_d(index));
                log_Phi(idx_node) = log_exp_x_plus_exp_y_vec(log_Phi(idx_node),
                                                             log_Phi(left_index) + log_Phi(right_index) + M_d(index));
            };
        }
        
        uvec idx_node = F(l);
        uvec compact_idx_node = compact_idx_map(idx_node);
        log_Phi(compact_idx_node) -= log(num_child_double(idx_node));
        
        log_Phi(compact_idx_node) = log_exp_x_plus_exp_y_vec(log(eta) + log_p0(idx_node), log(1 - eta) + log_Phi(compact_idx_node));
        
    };
    
    double a = log_Phi(compact_idx_map(0));
    return a;
}

mat class_tree::hyper2par(const vec &x)
{
    // x = [log.beta, log.C_2, qlogist(eta), log.alpha, log.tau, log.sigma]: length 6
    vec rho(total_level + 1), tau(total_level + 1), eta(total_level + 1);
    for (int i=0; i < total_level + 1; i++)
    {
        rho(i) = min(1.0, pow(2.0, -exp(x[0]) * i) * exp(x[1]));
        tau(i) = pow(2.0, -exp(x[3]) * i) * exp(x[4]);
        eta(i) = 1.0/(1.0 + exp(-x[2]));
    };
    
    eta(total_level) = exp(x[5]);
    
//     if (x.n_elem > 5)
//     {
//         eta(total_level) = exp(x[5]);
//     } else 
//     {
//         eta(total_level) = 1.0; // sigma = 1 (observation rescaled) 
//     }; 
    
    mat par(total_level + 1, 3);
    par.col(0) = rho; par.col(1) = tau; par.col(2) = eta;
    return par;
};


vec class_tree::fit_tree(const mat& par, const vec& shift_vec)
{
    
    // parameters: 3 rows by (total_level + 1): last element is sigma
    vec rho_vec = par.col(0);
    vec tau_vec = par.col(1);
    vec eta_vec = par.col(2);
    double sigma = eta_vec(total_level);
    eta_vec.shed_row(total_level);
    
    // generate phi map
    // N_prime: log_phi, log_p0
    // pair_N: M_d, rho_d, post_lambda_d
    
    vec log_rho = log(rho_vec);
    vec log_rho_c = log(1 - rho_vec);
    
    vec w_d(N, fill::zeros), ss_adjust(N, fill::zeros);
    
    if (all(shift_vec == zeros(shift_vec.n_elem)))
    {
        w_d = w_d_orig;
        ss_adjust = ss_adjust_orig;
    } else {
        vec shift_X = circshift_vectorize(obs, dimension, shift_vec);
        vec y(N, fill::zeros), ss(N, fill::zeros);
        y(F(total_level)) = shift_X;
        ss(F(total_level)) = shift_X % shift_X;
        
        for (int l = total_level - 1; l >= 0;l --){ // current level
            uvec index = flag_first(l + m * (total_level + 1)); // select the current node has dth cut & first cut
            uvec left_index = left_idx_d(index);
            uvec right_index = right_idx_d(index);
            uvec idx_node = node_idx_d(index);
            y(idx_node) = y(left_index) + y(right_index);
            ss(idx_node) = ss(left_index) + ss(right_index);
        }
        
        ss_adjust = ss - ((y % y) / n_A);
        w_d = y(left_idx_d) - y(right_idx_d);
        w_d.each_col() /= sqrt(n_A(node_idx_d));
    }
    
    vec log_Phi = get_phi_map(par, w_d, ss_adjust);
    
    uvec invert_compact_idx_map = find(num_child > 0);
    
    //    uword test = all(compact_idx_map(invert_compact_idx_map) == linspace(0, N_prime - 1, N_prime));
    vec log_p0 = (n_A(invert_compact_idx_map) - 1) * (log_inv_sqrt_2pi - log(sigma)) - (0.5 * (1/sigma) * (1/sigma)) * ss_adjust(invert_compact_idx_map);
    
    vec log_eta = log(eta_vec);
    vec log_eta_c = log(1 - eta_vec);
    
    
    vec post_eta(N_prime);
    post_eta = log_eta(level(invert_compact_idx_map)) + log_p0 - log_Phi;
    
    vec  post_lambda_d(pair_N, fill::zeros), M_d(pair_N), post_rho_d(pair_N);
    for (int i = 0; i < pair_N; i++)
    {
        double part1 = log_rho(pair_level(i)) + log_ratio_function(w_d(i) / sigma, tau_vec(pair_level(i)));
        double part2 = log_exp_x_plus_exp_y(part1, log_rho_c(pair_level(i)));
        M_d(i) = part2 + normal_logpdf(w_d(i), sigma);
        post_rho_d(i) = part1 - part2;
        
        if (num_child_double(node_idx_d(i)) > 1) // post_lambda_d: last level meaningness; 2nd last = 1 (log = 0); generally if it only has one child, the prob = 1 (log = 0)
        {
            post_lambda_d(i) = log_eta_c(pair_level(i)) - log(num_child_double(node_idx_d(i))) + M_d(i) +
            log_Phi(compact_idx_map(left_idx_d(i))) + log_Phi(compact_idx_map(right_idx_d(i))) - log_Phi(compact_idx_map(node_idx_d(i))) - log(1 - exp(post_eta(compact_idx_map(node_idx_d(i)))));
        };
        
    };
    
    vec psi_0(N, fill::zeros), varphi_0_kernel(N, fill::zeros), varphi(N, fill::zeros);
    psi_0(0) = 1 - exp(post_eta(compact_idx_map(0)));
    
    varphi_0_kernel(0) = rescale_sum;
    varphi(0) = rescale_sum ; // modify later on
    
    // use the raw scale (not "log") for these three maps:
    // for psi_0: probabilities (possible to use log scale but let's see how it works)
    // for varphi_0_kernel and varphi: could be negative so no way to use log(e^x + e^y)
    
    // last level is = 0
    for (int l = 0; l < total_level - 1; l ++) // l is the current level as parents: update (l + 1)th level
    {
        // recursive update: use flag_not_first to sum across dimension 'd'
        uvec index = flag_not_first(l + m * (total_level + 1)); // index in N_pair
        uvec node = node_idx_d(index);
        vec multiplier = 1 - exp(post_eta(compact_idx_map(node)));
        
        uvec node_left_child = left_idx_d(index);
        uvec node_right_child = right_idx_d(index);
        vec addup = psi_0(node) % exp(post_lambda_d(index));
        
        psi_0(node_left_child) += addup;
        psi_0(node_right_child) += addup;
        
        //one-time update: use "F(l)" since "flat_not_first" has replicates
        node = F(l + 1);
        multiplier = 1 - exp(post_eta(compact_idx_map(node)));
        psi_0(node) %= multiplier;
        
    };
    
    vec post_mean_map(pair_N);
    post_mean_map = sigma * (exp(post_rho_d) % mu_one(w_d / sigma, tau_vec(pair_level)));
    
    for (int l = 0; l < total_level; l ++) // l is the current level as parents: update (l + 1)th level
    {
        
        // recursive update: use flag_not_first to sum across dimension 'd'
        uvec index = flag_not_first(l + m * (total_level + 1)); // index in N_pair
        uvec node = node_idx_d(index);
        vec multiplier = 1 - exp(post_eta(compact_idx_map(node)));
        
        uvec node_left_child = left_idx_d(index);
        uvec node_right_child = right_idx_d(index);
        
        vec add_part1 = sqrt(0.5) * (exp(post_lambda_d(index)) % varphi_0_kernel(node) % multiplier);
        vec add_part2 = sqrt(0.5) * (exp(post_lambda_d(index)) % post_mean_map(index) % psi_0(node));
        varphi_0_kernel(node_left_child) += (add_part1 + add_part2);
        varphi_0_kernel(node_right_child) += (add_part1 - add_part2);
        
        vec add_part3 = sqrt(0.5) * (varphi(node) - varphi_0_kernel(node) % multiplier) / num_child_double(node);
        varphi(node_left_child) += (add_part1 + add_part2 + add_part3);
        varphi(node_right_child) += (add_part1 - add_part2 + add_part3);
        
    };
    
    vec BMA = varphi(F(total_level));
//    BMA.set_size( dimension );
    
    return BMA; 
//     //MAP tree
//     vec log_kappa(N_prime, fill::zeros); // update using the same flow of log_Phi
//     uvec idx_pair(N_prime); // add: if selected: the idx in pair_N selected
//     
//     // 2nd last level: direction - only one possibility; log_kappa = 0
//     int l = total_level - 1;
//     uvec index = flag_first(l + m * (total_level + 1));
//     uvec idx_node = compact_idx_map(node_idx_d(index));
//     idx_pair(idx_node) = index;
//     
//     for (int l = total_level - 2; l >=0; l--)
//     {
//         for (int d = 0; d < m; d ++)
//         {
//             uvec index = flag_first(l + d * (total_level + 1)); // if the current node has dth cut & first cut
//             uvec left_index = compact_idx_map(left_idx_d(index));
//             uvec right_index = compact_idx_map(right_idx_d(index));
//             uvec idx_node = compact_idx_map(node_idx_d(index));
//             log_kappa(idx_node) = (log_kappa(left_index) + log_kappa(right_index) + post_lambda_d(index));
//             idx_pair(idx_node) = index;
//             
//             index = flag_not_first(l + d * (total_level + 1)); // if the current node has more than dth cut
//             left_index = compact_idx_map(left_idx_d(index));
//             right_index = compact_idx_map(right_idx_d(index));
//             idx_node = compact_idx_map(node_idx_d(index));
//             vec candidate = log_kappa(left_index) + log_kappa(right_index) + post_lambda_d(index);
//             uvec idx = find(log_kappa(idx_node) < candidate);
//             log_kappa(idx_node(idx)) = candidate(idx);
//             idx_pair(idx_node(idx)) = index(idx);
//             
//         };
//         
//     };
//     
//     //this block is potentially useful for visulization
//     //    // if selected: left child and right child node index
//     //    uvec left_node = left_idx_d(idx_pair);
//     //    uvec right_node = right_idx_d(idx_pair);
//     //    // for each node (has children): if selected which direction
//     //    uvec direction = pair_direction(idx_pair);
//     
//     
//     // a tree: has 2^(total_level + 1) - 1 elements; same as R (last 2^total_level elements meaningsless)
//     // use field to be consistent with R coding
//     // (which pair selected, R of the selected pair/node, z map)
//     field<uvec> pair_in(total_level), node_R(total_level);
//     // uvec indicator(N_prime, fill::ones), R(N_prime); // pruned or not: generate at the last step
//     pair_in(0) = idx_pair(0);
//     node_R(0) = (post_eta(compact_idx_map(node_idx_d(pair_in(0)))) > log(0.5));
//     
//     for (int l = 0; l < total_level - 1; l ++) // current level: decide best 'd' and children nodes included
//     {
//         uvec idx = pair_in(l);
//         
//         uvec left_pair_in = idx_pair(compact_idx_map(left_idx_d(idx)));
//         uvec right_pair_in = idx_pair(compact_idx_map(right_idx_d(idx)));
//         pair_in(l + 1) = merge_left_right(left_pair_in, right_pair_in);
//         
//         uvec multiplier_temp = merge_left_right(node_R(l), node_R(l));
//         node_R(l + 1) = ((post_eta(compact_idx_map(node_idx_d(pair_in(l + 1)))) > log(0.5)) || multiplier_temp);
//         
//     };
//     
//     // calculate z map: has total_level's elements
//     field<vec> z_map(total_level);
//     for (int l = 0; l < total_level; l ++)
//     {
//         z_map(l) = post_mean_map(pair_in(l)) % (1 - node_R(l));
//     };
//     
//     // calcualte MAP estimates
//     vec a(1);
//     a(0) = rescale_sum;
//     for (int l = 0; l < total_level; l ++)
//     {
//         vec left = sqrt(0.5) * (a + z_map(l)); vec right = sqrt(0.5) * (a - z_map(l));
//         a = merge_left_right(left, right);
//     }
//     
//     uvec left_node = left_idx_d(pair_in(total_level - 1));
//     uvec right_node = right_idx_d(pair_in(total_level - 1));
//     uvec node_last = merge_left_right(left_node, right_node);
//     
//     vec MAP_temp(N);
//     MAP_temp(node_last) = a;
//     mat MAP = MAP_temp(F(total_level));
//     MAP.set_size( size(obs) );
//     
//     field<mat> ret(2);
//     ret(0) = BMA; ret(1) = MAP;
//     return ret;
    
}



vec class_tree::fit_tree_cs(const mat& par, int step)
{
    
    vec global_BMA(obs.n_elem, fill::zeros);
   // global_BMA.set_size( size(obs) );
   // global_MAP(size(obs), fill::zeros); global_MAP.set_size( size(obs) );

    // two-dimensional coding:
//    vec a1 = regspace<vec>(-step, step);
//    mat a2 = repmat(a1, 1, 2*step + 1);
//    mat a3 = repmat(a1.t(), 2*step + 1, 1);
//    mat all_shift(2, (2*step + 1) * (2 * step + 1));
//    all_shift.row(0) = vectorise(a2).t();
//    all_shift.row(1) = vectorise(a3).t();
//    
// 
    
    // int count = 0;
    uword m = dimension.n_elem;
    uword num_shift = pow((2 * step + 1), m);
    uvec shift_basis(m);
    shift_basis.fill(2 * step + 1);
    vec center_basis(m);
    center_basis.fill(step);
    
    for (uword ith_shift = 0; ith_shift < num_shift; ith_shift ++)
    {
        vec shift_vec = conv_to< vec >::from(ind2sub(ith_shift, shift_basis)) - center_basis;        
        vec est = fit_tree(par, shift_vec);
        global_BMA += (circshift_vectorize(est, dimension, -shift_vec) / num_shift);
        // global_MAP += (circshift(est(1), -all_shift.col(ith_shift)) / num_shift);
        // count ++;
    }
    
    
    return global_BMA; 
//    global_BMA /= count;
//    global_MAP /= count;
    
//     field<mat> ret(2);
//     ret(0) = global_BMA; ret(1) = global_MAP;
//     return ret;
    
}

vec class_tree::get_M_d_map(const mat& par)
{
    
    // parameters: 3 rows by (total_level + 1): last element is sigma
    vec rho_vec = par.col(0);
    vec tau_vec = par.col(1);
    vec eta_vec = par.col(2);
    double sigma = eta_vec(total_level);
    eta_vec.shed_row(total_level);
    
    vec M_d(pair_N);
    vec log_rho = log(rho_vec);
    vec log_rho_c = log(1 - rho_vec);
    
    vec w_d(N, fill::zeros), ss_adjust(N, fill::zeros);
    w_d = w_d_orig;
    ss_adjust = ss_adjust_orig;
    
    for (int i = 0; i < pair_N; i++)
    {
        M_d(i) = log_exp_x_plus_exp_y(log_rho(pair_level(i)) + log_ratio_function(w_d(i) / sigma, tau_vec(pair_level(i))),
                                      log_rho_c(pair_level(i)) ) + normal_logpdf(w_d(i), sigma);
        
    };
    
    return M_d; 
}

vec class_tree::get_lambda_map(const mat& par)
{
    
    // parameters: 3 rows by (total_level + 1): last element is sigma
    vec rho_vec = par.col(0);
    vec tau_vec = par.col(1);
    vec eta_vec = par.col(2);
    double sigma = eta_vec(total_level);
    eta_vec.shed_row(total_level);
    
    // generate phi map
    // N_prime: log_phi, log_p0
    // pair_N: M_d, rho_d, post_lambda_d
    
    vec log_rho = log(rho_vec);
    vec log_rho_c = log(1 - rho_vec);
    
    vec w_d(N, fill::zeros), ss_adjust(N, fill::zeros);
    // in 'fit_tree' function, we allow cycle spinning to obtain various 'w_d' and 'ss_adjust'; 
        w_d = w_d_orig;
        ss_adjust = ss_adjust_orig;
    
    vec log_Phi = get_phi_map(par, w_d, ss_adjust);
    
    uvec invert_compact_idx_map = find(num_child > 0);
    
    //    uword test = all(compact_idx_map(invert_compact_idx_map) == linspace(0, N_prime - 1, N_prime));
    vec log_p0 = (n_A(invert_compact_idx_map) - 1) * (log_inv_sqrt_2pi - log(sigma)) - (0.5 * (1/sigma) * (1/sigma)) * ss_adjust(invert_compact_idx_map);
    
    vec log_eta = log(eta_vec);
    vec log_eta_c = log(1 - eta_vec);
    
    
    vec post_eta(N_prime);
    post_eta = log_eta(level(invert_compact_idx_map)) + log_p0 - log_Phi;
    
    vec  post_lambda_d(pair_N, fill::zeros), M_d(pair_N), post_rho_d(pair_N);
    for (int i = 0; i < pair_N; i++)
    {
        double part1 = log_rho(pair_level(i)) + log_ratio_function(w_d(i) / sigma, tau_vec(pair_level(i)));
        double part2 = log_exp_x_plus_exp_y(part1, log_rho_c(pair_level(i)));
        M_d(i) = part2 + normal_logpdf(w_d(i), sigma);
        post_rho_d(i) = part1 - part2;
        
        if (num_child_double(node_idx_d(i)) > 1) // post_lambda_d: last level meaningness; 2nd last = 1 (log = 0); generally if it only has one child, the prob = 1 (log = 0)
        {
            post_lambda_d(i) = log_eta_c(pair_level(i)) - log(num_child_double(node_idx_d(i))) + M_d(i) +
            log_Phi(compact_idx_map(left_idx_d(i))) + log_Phi(compact_idx_map(right_idx_d(i))) - log_Phi(compact_idx_map(node_idx_d(i))) - log(1 - exp(post_eta(compact_idx_map(node_idx_d(i)))));
        };
        
    };
    
    return post_lambda_d; 
}




