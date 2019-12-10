% SparseUpdateOneLevel updates (mother) wavelet coefficients if the next
% level father wavelet coefficients change sparsely
% father coefficients 'father_child' are updated using (idx_child,
% value_child), i.e., the new father coef of the child level become
% father_child + value_child at idx_child
% (approx_coef, detail_coef)
% 'level' is the updating level; number of wavelet coefs = 2^level
function [p_idx, p_diff] = SparseUpdateOneLevel(c_idx, c_diff, operator, level)
% idx_child starts with length length(g)/2, affecting at most length(g) - 1

% all vectors are column vectors

if c_idx < 0
    error('c_idx must be nonnegative');
end
% c_idx takes value from 0, 1, 2, ..., (2^c_level - 1)


    c_idx_last = mod(c_idx + length(c_diff) - 1, 2^level);
    flag = [mod(c_idx, 2), mod(c_idx_last, 2)];
    
    % expand c_diff into size-2 blocks
    c_bin_diff = c_diff;
    
    if (flag(1) == 1)
        c_bin_diff = [0; c_bin_diff];
    end
    
    if (flag(2) == 0)
        c_bin_diff = [c_bin_diff; 0];
    end
    
    %% affects to the parent level
    % block_idx in the children levels is affected
    block_idx = floor(c_idx/2);
    l = length(operator)/2; % operator must have even length
    p_idx = mod(block_idx - l + 1, 2^(level - 1));
    
    b = length(c_bin_diff)/2;
    p_length = b + l - 1;
    p_diff = zeros([p_length, 1]);
    
    % utilize the fact that the number of affected blocks b is at most l
    for i = 1:(b + l - 1) % number of overlaped items
        for j = max(1, i - l + 1):min(i, b)
            % each overlapyed item uses jth block * (l + j - i)th operator block
            p_diff(i) = p_diff(i) + ...
                + operator(2 * (l + j - i) - 1) * c_bin_diff(2*j - 1) + ...
                + operator(2 * (l + j - i)) * c_bin_diff(2*j);
        end
    end

end 




