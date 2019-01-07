# WARP
Wavelets with adaptive recursive partitioning applied to image reconstruction (both 2D and 3D). Matlab users can directly use the code for real data analysis. The C++ source code is also available, which is portable to`R` via `RcppArmadillo`. 

# Instructions
- Install the c++ library `armadillo`(http://arma.sourceforge.net) to your computer; 
- Let `path` be the directory containing the header file of `armadillo`; 
- Compile the c++ code for matlab by running the code `mex_this(path)`. 

# 2D example

The following code is available in the file `Demo.m`. It usually takes 20sec if ran in a Macbook Pro with usual specifications, which increases by 1min if the step size in cycling spinning (`step`) is changed to 5.  

```

rng(0);
load('lena_data.mat');
obs_true = double(lena) ./ 255; % ground truth
MSE = @(x) mean(mean((x - obs_true).^2, 1));

sigma = 0.2;
obs_raw = obs_true + randn(size(obs_true, 1)) .* sigma; % 2D observation
obs = obs_raw(:); % vectorize obs
dimension = size(obs_raw)'; % obtain size information as a column vector

tic
% default method to select parameters in the model
hyper0 = hyper_default(obs(:), dimension);
t1 = toc;
% Bayesian model averaging without cycle spinning
BMA_no_cs = treeFit(obs, dimension, hyper0, 0);
BMA_no_cs = reshape(BMA_no_cs, dimension'); 
t2 = toc;
% BMA with cycle spinning;
step = 1; % increase `step` will lead to better performance (smaller MSEs) but takes longer time. 
BMA_cs = treeFit(obs, dimension, hyper0, step);
BMA_cs = reshape(BMA_cs, dimension'); 
t3 = toc;

figure; 
subplot(2, 2, 1)
imshow(obs_true); title('True'); 
subplot(2, 2, 2)
imshow(obs_raw); 
title(sprintf('Noisy observation \n MSE = %.4f', MSE(obs_raw))); 
subplot(2, 2, 3)
imshow(BMA_no_cs); 
title(sprintf('WARPed Haar without CS \n MSE = %.4f | Time = %.1fs', MSE(BMA_no_cs), t2)); 
subplot(2, 2, 4)
imshow(BMA_cs); 
title(sprintf('WARPed Haar with %d shifts \n MSE = %.4f | Time = %.1fs', (2 * step + 1)^2, MSE(BMA_cs), t3 - t2 + t1));
```

Although we do not recommend it, if you'd like to use a full optimization to select hyperparameters rather than the default method implemented in `hyper_default`, use the following code to obtain `hyper0`.

```
idx = (1:2:numel(obs));
a = (obs(idx) - obs(idx + 1)) ./ sqrt(2);
sigma_hat = mad(a(:), 1) * 1.4826;
x0 = [0, 0, 0, 0, 0]';
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
f = @(x) -treeLikelihood(obs,[x; log(sigma_hat)]);
[x, fval] = fminunc(f,x0, options);
hyper0 = [x; log(sigma_hat)];
```
