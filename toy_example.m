%%%% toy example  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% Global settings %%%%%%
str = computer;
% use 'mex_this.m' to compile the cpp code for matlab if neither works
addpath('source_3D_matlab_mac');
addpath('source_SMC');
% addpath('WARP');
rng(0);

%% Generate the toy stripe flag %%%%%

mat = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9;
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8;
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9;
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8];
   sigma=0.01; mat = mat + randn(size(mat))*sigma;
%imagesc(mat) %show your original image
dimension = size(mat);
n = numel(mat);

%% CARP for compression %%%%%
sigma_hat = 0.01;% the adjustable parameter related to CR
n_tree = 2;% number of random trees
rand_seed = 2019; %random seed used by armadillo in the function "DrawPosition"
% consider to change rand_seed if run multiple replications
hyper1 = hyper_3D_default(mat(:), dimension', true, sigma_hat);
[MAP_0, MAP_fit] = MAP(mat(:), dimension', hyper1);
% layout of 'smp_all':
direction_1 = MAP_0(:, 1);
pruning_1 = MAP_0(:, 2);
% shorten direction
pruning_opt = direction_1;
pruning_opt(pruning_1==1) = 3;
% find all the children
J = log2(prod(dimension));
level = J+1;
L = length(pruning_opt);
for l=1:(level-2)
    node_l = pruning_opt(2^(l-1):(2^l-1));
    indx = find(node_l==3);
    if isempty(indx)==0
        d = length(indx);
        for t=1:d
            node_idx = 2^(l-1)-1+indx(t);
            %pruning_0(node_idx) = 4;
            child_indx = [];
            for j=1:(level-l-1)
                child_indx = [child_indx, (node_idx*2^j):(node_idx*2^j+2^j-1)];
            end
            pruning_opt(child_indx) = 4;
        end
    end
end
ShortenPruning = double(pruning_opt(pruning_opt ~= 4));
LengthRateShort = length(ShortenPruning)/L;
ShortenDirection = ShortenPruning;
ShortenDirection(ShortenDirection==3) = -1;
prun_num = sum(ShortenDirection==-1);
[position_new, prun_block] = dir2pos_short(ShortenDirection, prun_num, dimension);
rec_mat = mat;
comp_mat = mat(:);
total_block = [];
total_block_avg = zeros(prun_num,1);
for t=1:prun_num
    block_t = prun_block{t,1};
    total_block = [total_block;reshape(block_t,[numel(block_t),1])];
    total_block_avg(t) = mean(rec_mat(block_t(:)));
    rec_mat(block_t(:)) = total_block_avg(t);
end
comp_mat(total_block) = [];
comp_mat = [comp_mat;total_block_avg];
comp_mat_size = length(comp_mat);
%% wavelet decompostion & reconstruction
%rec_mat_vec = rec_mat(position_new);
[rec_c0, rec_wave_level] = wavedec(comp_mat, log2(n), 'db1');
% soft threshold
r1=0.0001;
coeff_mat_rec = rec_c0;
coeff_mat_rec(abs(rec_c0)<r1) = 0;
coeff_mat_rec(rec_c0>r1) = coeff_mat_rec(rec_c0>r1)-r1;
coeff_mat_rec(rec_c0<-r1) = coeff_mat_rec(rec_c0<-r1)+r1;
%ytsoft = wthresh(rec_c0,'s',r1);
rec_mat_vec_final = waverec(coeff_mat_rec,rec_wave_level,'db1');
final_mat = rec_mat;
block_num = length(prun_block);
for kk =1:block_num
    final_mat(prun_block{kk,1}) =...
        rec_mat_vec_final(comp_mat_size-block_num+kk);
end
rec = final_mat;

%% Huffman encoding
coeff_mat_rec_ed = coeff_mat_rec;
coeff_mat_rec_ed(coeff_mat_rec_ed==0) = [];
coeff_mat_rec_ed = floor((coeff_mat_rec_ed-min(coeff_mat_rec_ed))/...
    (max(coeff_mat_rec_ed)-min(coeff_mat_rec_ed))*255);
symbols = unique(coeff_mat_rec_ed);
Cat = categorical(coeff_mat_rec_ed,symbols);
prob = histcounts(Cat)/length(coeff_mat_rec_ed);
prob = prob/(sum(prob));
dict = huffmandict(symbols, prob);
enco = huffmanenco(coeff_mat_rec_ed, dict);
FinalCompressedImage = numel(de2bi(enco));
% Calculate the compression rate
binarySig = de2bi(floor(255*mat));
seqLen = numel(binarySig);
compressRate  = seqLen/(FinalCompressedImage+log2(3)*comp_mat_size);
%compres = 8*n/(8*6+log2(3)*comp_mat_size);

%% Calculate the PSNR and MSE
D = abs(double(mat)-double(final_mat)).^2;
mse  = sum(D(:))/prod(dimension);
psnr = 10*log10(1^2/mse);
fprintf('Compression ratio is %6.3f.\n', compressRate)
fprintf('PSNR value is %6.3f.\n', psnr)
imagesc(final_mat)