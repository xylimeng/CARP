function [output, direction] = DrawChildrenNodes(rank_node, t, lambda_mat, rank_left_child, rank_right_child, keep_all)
% draw t-step children nodes of rank_node (scalar)
% 'lambda_mat', 'rank_left_child', 'rank_right_child' are global variables

if nargin < 6
    keep_all = 1; % keep all levels or just the very last child level 
end 
N = size(lambda_mat, 1);
output = rank_node; 
direction = null(1); 
for ith_step = 1:t
    divide_prob = lambda_mat(rank_node, :);
    flag = (rand([size(divide_prob, 1), 1]) < divide_prob(:, 1)); % if 1, then diretion 1; if 0, then direction 2
    which_direction = 2 - flag;
    % direction_in{l + 1} = which_direction; % row: rank_node; col: which_direction
    idx_selected = (uint64(which_direction) - 1) * N + rank_node;
    
    rank_node_next_level = cat(1, rank_left_child(idx_selected)', rank_right_child(idx_selected)'); % rbind - 2 rows
    rank_node = rank_node_next_level(:);
    output = [output; rank_node]; 
    direction = [direction; which_direction]; 
end
if keep_all == 0 
    output = rank_node; 
end 
end