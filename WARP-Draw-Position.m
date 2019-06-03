%% Compare JPEG & JPEG2000 with WARP in 3D image
%% fix a compression rate, compare the PSNR
%% Rongjie Liu @ 03/2019 @ Rice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WARP for reducing dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global settings
str = computer;
% either mac or linux
% use 'mex_this.m' to compile the cpp code for matlab if neither works
if strcmp(str, 'MACI64')
    addpath('source_3D_matlab_mac');
else
    addpath('source_3D_matlab_linux');
end
addpath('source_SMC');
addpath('WARP');
rng(0);
% R: the maximum fluctuation in the input image data type
% R = 1; % double-precision floating-point data type
% R = 255; % 8-bit unsigned integer data type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('brain.mat');
%X = imresize3(V,[64, 64, 64]);
%obs = double(floor(255*(X-min(X(:)))/(max(X(:))-min(X(:)))));
X = reshape(1:8, [2,2,2])/8;
obs = X;
dimension = size(obs);
hyper0_t = hyper_3D_default(obs(:), dimension', false);

% number of random trees
n_tree = 2;
rand_seed = 2019; %random seed used by armadillo in the function "DrawPosition"
smp_all = DrawPosition(obs(:), dimension', hyper0_t, n_tree, ...
                            rand_seed);
% consider to change rand_seed if run multiple replications 

%% layout of 'smp_all': 
direction_all = smp_all(1:(numel(obs) - 1), :); 
pruning_all = smp_all(numel(obs): (2 * numel(obs) - 2), :); 
position_all = smp_all((2 * numel(obs) - 1):end, :); 

% summarize one tree
ith_tree = 1;
position = position_all(:, ith_tree);
oneD = obs(position + 1); % position returned by c++ starts from 0;
                          % add 1 for matlab 
