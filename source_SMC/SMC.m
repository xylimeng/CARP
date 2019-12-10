function [final_est, final_weight, est_all] = SMC(obs, dimension, num_particle, ESS, ...
    lambda_map, rank_left_child, rank_right_child, dividible, par, y, length_filter)


% global settings
total_level = log2(numel(obs));
num_node = 2 * dimension - 1;
N = prod(num_node);
sigma_hat = par(total_level + 1, 3);
dwtmode('per'); %periodization:this boundary mode makes c the same length as obs

% field 'n' stores various length information
% clear n;
n.obs = numel(obs);
% wavelet filter length = 2 * n.filter;
n.filter = length_filter; % 'db2' filter ('db1' is Haar)
nm = sprintf('db%d', n.filter);

[Lo, Hi, ~, ~] = wfilters(nm);
n.tree = num_particle; % number of trees
n.dict_node = N; % number of nodes in the dictionary
n.tree_size = n.obs - 1; % size of each tree = num of observations less 1
n.tree_level = total_level; % max level in each tree

% in cpp: code this efficiently
lambda_mat = zeros([n.dict_node, 2]);
lambda_mat((dividible > 0)) = exp(lambda_map);
% M_d_mat = zeros([n.dict_node, 2]);
% M_d_mat((dividible > 0)) = M_d_map;
log_ratio_lambda_mat = zeros([n.dict_node, 2]); % log(prior / posterior)
log_ratio_lambda_mat((dividible > 0)) = log(1/2) - lambda_map;

wave_coef = zeros([n.obs, n.tree]);
wave_coef(1,:) = y(1)/sqrt(n.obs);
% available node for branching after accounting for near atomic & pruning:
wave_open = ones(size(wave_coef)); % level 0, 1, ..., n.tree_level - 1 INCLUDING 1 father coef - same index as 'wave_coef'

% % 1 - available; 0 - taken
% node_available = ones([n.obs - 1, n.tree]);
% % starting idx and ending idx at each level in Matlab (starting from 0)
% node_available_start = 2.^(0:(n.tree_level - 1));
% node_available_end = 2.^(1:(n.tree_level)) - 1;
% % levels 0 to (level_cutoff) are taken
% node_available(1:node_available_end(level_cutoff + 1), :) = 0;


wave_keeper = repelem(0:(n.tree_level - 1), 2.^(0:(n.tree_level - 1))); % house keeper for wave_coef less father coef
wave_coef_logL = zeros([n.obs, n.tree]);

logL_target = zeros([n.tree, 1]);
diff_log_ratio_lambda = zeros([n.tree, 1]); % store difference of log_ratio_lambda at each move (dynamic)
diff_logL_target = zeros([n.tree, 1]); % = new - old
% diff_log_weight = zeros([n.tree, 1]); % increment of log(w) = log (ratio_lambda * ratio_likelihood)


% initialize
node_in = ones([1, n.tree]); % the root node
log_weight = zeros([n.tree, 1]);


% convinient utility functions
conv_DrawChildren = @(rank, step) DrawChildrenNodes(rank, ...
    step, lambda_mat, rank_left_child, rank_right_child);
% conv_DrawChildren_neat = @(rank, step) DrawChildrenNodes(rank, ...
%     step, lambda_mat, rank_left_child, rank_right_child, 0); % only keep the last level



%% starting from level 0: draw multiple steps simultaneously
% if level_cutoff < min(dimension(1), dimension(2)) then no need to
% consider nearly atomic nodes but possible has prunning
level_cutoff = length_filter;
keeper = repelem(0:level_cutoff, 2.^(0:level_cutoff))';
idx_cutoff = find(keeper == level_cutoff);
node_in_new = zeros([2^level_cutoff, n.tree]);

for ith_tree = 1:n.tree
    rank_node = node_in(1, ith_tree);
    [child_node, d] = conv_DrawChildren(rank_node, level_cutoff);
    % generate the likelihood ratio
    % temp_obs = sum of obs at each node / sqrt(number of pixels at each node)
    temp_obs = y(child_node(idx_cutoff)) ./ sqrt(2^(n.tree_level - level_cutoff));
    node_in_new(:, ith_tree) = child_node(idx_cutoff);
    
    % wavelet coefficients
    [c2, l2] = wavedec(temp_obs, level_cutoff, nm);
    
    for l = 0:(level_cutoff - 1)
        wave_coef((2^l:(2^(l + 1) - 1)) + 1, ith_tree) = detcoef(c2, l2, level_cutoff - l);
        wave_open((2^l:(2^(l + 1) - 1)) + 1, ith_tree) = wave_open((2^l:(2^(l + 1) - 1)) + 1, ith_tree) + 1;
    end
    
    for l = 1:(2^level_cutoff - 1)
        % logL_haar(ith_tree) = logL_haar(ith_tree) + M_d_mat(child_node(l), d(l));
        diff_log_ratio_lambda(ith_tree) = diff_log_ratio_lambda(ith_tree) + log_ratio_lambda_mat(child_node(l), d(l));
        wave_coef_logL(l + 1, ith_tree) = wavecoef2logL(wave_coef(l + 1, ith_tree), wave_keeper(l), par, sigma_hat);
        logL_target(ith_tree) = logL_target(ith_tree) + wave_coef_logL(l + 1, ith_tree);
    end
    log_weight(ith_tree) = diff_log_ratio_lambda(ith_tree) + logL_target(ith_tree);
end

weight = exp(log_weight - max(log_weight)); W = sum(weight(:));
weight_norm = weight ./ W;
log_W = log(W) + max(log_weight);

fprintf('ESS = %.2f \n', 1 / sum(weight_norm .^2));
if (1 / sum(weight_norm .^2) < ESS) % resample
    fprintf('resample needed | \n');
    flag_resample = zeros([n.tree, 1]);
    for ith_tree = 1:n.tree
        flag_resample(ith_tree) = find(mnrnd(1, weight_norm'));
    end
    log_weight = zeros([n.tree, 1]) + log_W - log(n.tree);
    
    % re-assign according to flag_resample;
    node_in_new = re_assign(node_in_new, flag_resample);
    wave_coef = re_assign(wave_coef, flag_resample);
    wave_coef_logL = re_assign(wave_coef_logL, flag_resample);
    diff_log_ratio_lambda = re_assign(diff_log_ratio_lambda, flag_resample);
    logL_target = re_assign(logL_target, flag_resample);
end

% fprintf(sprintf('ends: took %.2fs \n', toc));

%% next levels: node by node
% legend:
%   j - branching level (previous branching_level)
%   k - branching node (from 0 to 2^j - 1)

for j = level_cutoff:(n.tree_level - 1)
    tic
    fprintf(sprintf('%dth level starts: resample at ', j));
    child_level = j + 1;
    
    % pre-allocate level storage
    node_in = node_in_new;
    node_in_new = zeros([2^child_level, n.tree]);
    dummy_obs = y(node_in) ./ sqrt(2^(n.tree_level - j));  %obs at jth level - old
    update_obs = zeros([2^j, n.tree]); % dynamically changed obs at jth level
    
    % determine nearly-atomic nodes
    
    for k = 0:(2^j - 1)
        for ith_tree = 1:n.tree
            rank_node = node_in(k + 1, ith_tree);
            [child_node, d] = conv_DrawChildren(rank_node, 1);
            temp_obs = y(child_node(2:3)) ./ sqrt(2^(n.tree_level - child_level));
            node_in_new((2*k + 1):(2*k + 2), ith_tree) = child_node(2:3);
            diff_log_ratio_lambda(ith_tree) = log_ratio_lambda_mat(rank_node, d);
            [a, d, c_idx] = dwt_sparse(j + 1, 2*k, temp_obs, Lo, Hi);
            
            % changes in the branching level 'approximation coef'
            idx = (mod(c_idx:(c_idx + length(d) - 1), 2^j))';
            update_obs(idx + 1) = d + update_obs(idx + 1);
            c_diff = update_obs(idx + 1) - dummy_obs(idx + 1);
            dummy_obs(idx + 1) = update_obs(idx + 1);
            
            % update wavelet coef at all previous levels
            [output_idx, output_diff, output_keeper] = SparseWavelet(j, c_idx, c_diff, Lo, Hi);
            % combine wavelet coef changes
            output_idx = [output_idx; idx + 2^j];
            output_diff = [output_diff; a];
            output_keeper = [output_keeper; zeros([length(idx), 1]) + j];
            % update likelihood
            ss = 0;
            for ii = 1:length(output_idx)
                wave_coef(output_idx(ii) + 1, ith_tree) = wave_coef(output_idx(ii) + 1, ith_tree) + output_diff(ii);
                new_logL = wavecoef2logL(wave_coef(output_idx(ii) + 1, ith_tree), output_keeper(ii), par, sigma_hat);
                ss = ss + new_logL - wave_coef_logL(output_idx(ii) + 1, ith_tree);
                wave_coef_logL(output_idx(ii) + 1, ith_tree) = new_logL;
            end
            diff_logL_target(ith_tree) = ss;
            
        end
        
        % update log_weight
        log_weight = log_weight + diff_logL_target + diff_log_ratio_lambda;
        
        weight = exp(log_weight - max(log_weight));
        W = sum(weight(:)); log_W = log(W) + max(log_weight);
        
        weight_norm = weight ./ W; 
        if (1 / sum(weight_norm .^2) < ESS) % resample
            fprintf('%d | ', k);
            % fprintf('ESS = %.2f : resample needed | \n', 1 / sum(weight_norm .^2));
            flag_resample = zeros([n.tree, 1]);
            for ith_tree = 1:n.tree
                flag_resample(ith_tree) = find(mnrnd(1, weight_norm'));
            end
            log_weight = zeros([n.tree, 1]) + log_W - log(n.tree);
            % re-assign according to flag_resample;
            node_in = re_assign(node_in, flag_resample);
            node_in_new = re_assign(node_in_new, flag_resample);
            wave_coef = re_assign(wave_coef, flag_resample);
            wave_coef_logL = re_assign(wave_coef_logL, flag_resample);
            logL_target = re_assign(logL_target, flag_resample);
        end
    end
    fprintf(sprintf('\n     ends: took %.2fs \n', toc));
end


%% Reconstruct images
est_all = zeros([numel(obs), n.tree]);

tic
for ith_tree = 1:n.tree
    [i, j] = ind2sub(num_node, node_in_new(:, ith_tree));
    position = sub2ind(dimension, i - dimension(1) + 1, j - dimension(2) + 1);
    obs_vec = obs(position);
    % wavelet denoising: use functions 'wavedec' and 'waverec'
    [c, wave_level] = wavedec(obs_vec, log2(length(obs_vec)), nm);
    
    %% calculate the likihoods and their ratios
    % tau - should have length (total_level) NOT total_level + 1
    % last element of tau is not useful - double check
    
    tau_level = repelem(par((1:n.tree_level), 2), wave_level(2:(n.tree_level + 1)));
    rho_level = repelem(par((1:n.tree_level), 1), wave_level(2:(n.tree_level + 1)));
    
    w_d = c(2:end); % wavelet coefficients
    nlogL = zeros([total_level, 2]);
    for l = 0:(total_level - 1)
        nlogL1 = -log(1 - par(l + 1, 1)) + normlike([0, sigma_hat],w_d(2^l:(2^(l+1) - 1)));
        nlogL2 = -log(par(l + 1, 1)) + normlike([0, sqrt(1 + par(l + 1, 2)) * sigma_hat], w_d(2^l:(2^(l+1) - 1)));
        nlogL(l+1,:) = [nlogL1, nlogL2];
    end
    logL = -nlogL;
    
    log_weight(ith_tree) = log_weight(ith_tree) + sum(log_expey(logL(:,1), logL(:,2))) - logL_target(ith_tree);
    
    %% smoothing - target wavelet - NOT haar wavelet
    % rho_tilde - posterior spike probability
    % = \frac{1}{ 1 + (1/rho - 1) * M_d_0 /M_d_1 }
    M_d_ratio = sqrt(1 + tau_level) .* exp(- w_d.^2 ./ (2*sigma_hat^2) ./ (1 + 1./tau_level));
    rho_post = 1 ./ (1 + (1./rho_level - 1) .* M_d_ratio);
    
    post_mean = rho_post .* w_d ./ (1 + 1./tau_level);
    est_target = waverec([c(1); post_mean]', wave_level, nm);
    
    est = zeros(dimension);
    est(position) = est_target;
    
    est_all(:, ith_tree) = est(:);
end

final_weight = normalize(log_weight);
final_est = sum(bsxfun(@times, est_all, final_weight), 2);

end







