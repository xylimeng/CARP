function [idx, c, l] = SparseWavelet(level, c_idx, c_diff, Lo, Hi)
% SPARSEWAVELET obtains wavelet coefficients for sparse observations; stored sparsely
% INPUT:
% length of observation = 2^level;
% a wavelet vector = (father coef, mother coef's);
% c_idx is the index of the first observation that is changed (on the observation domain)

c = []; idx = []; l = [];

x = c_diff(:)'; % row vector
cutoff_level = ceil(log2( 2 * length(Lo))); % such that 2 * length(filter) <= 2^cutoff_level


for i = (level - 1):-1:(cutoff_level - 1) % ith level to be updated
    [x, d, c_idx] = dwt_sparse(i + 1, c_idx, x, Lo, Hi);
    % idx in the original wavelet vector: father coef is included
    idx_update = mod(c_idx:(c_idx + length(x) - 1), 2^i) + 2^i;
    c = [d c]; idx = [idx_update idx]; l = [repelem(i, length(x)) l];
end

% sort x : fully update at coarse levels
i = cutoff_level - 1;
order_x = zeros([1, 2^i]); % row vector
order_idx = mod((c_idx:(c_idx + length(x) - 1)), 2^i);
order_x(order_idx + 1) = x;
% the global dwtmode is 'per'
[d, keeper] = wavedec(order_x, cutoff_level - 1, Lo, Hi);
c = [d c]; idx = [0:(2^i - 1) idx]; l = [0 repelem(0:(i - 1), keeper(2:(i + 1))) l];

% column vector output
c = c'; idx = idx'; l = l';
end
