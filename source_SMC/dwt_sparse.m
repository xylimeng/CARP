% dwt for sparse observations to mimic matlab 'dwt' function
% SparseUpdateOneLevel updates (mother) wavelet coefficients if the next
% level father wavelet coefficients change sparsely
% father coefficients 'father_child' are updated using (idx_child,
% value_child), i.e., the new father coef of the child level become
% father_child + value_child at idx_child
% (approx_coef, detail_coef)
% 'level' is the updating level; number of wavelet coefs = 2^level
function [a, d, c_idx] = dwt_sparse(level, c_idx, c_diff, Lo, Hi)
% idx_child starts with length length(g)/2, affecting at most length(g) - 1
% c_idx takes value from 0, 1, 2, ..., (2^c_level - 1)

% vectors could be both col or row vectors 
x = c_diff; 
% Compute sizes and shape.
lf = length(Lo);
lx = length(x);
lenEXT = lf/2; 

% Figure out: 'first', 'last', 'idx' and 'c_idx' 
first = mod(lenEXT + c_idx, 2) + 1; %matlab idx  -> cpp idx = first - 1
len_full_conv = lx + lf - 1; % matlab last idx 
c_idx = mod(ceil((c_idx - lenEXT)/2), 2^(level - 1)); 

% Compute coefficients of approximation.
z = wconv1(x,Lo,'full');
a = z(first:2:len_full_conv);

% Compute coefficients of detail.
z = wconv1(x,Hi,'full');
d = z(first:2:len_full_conv);

end
