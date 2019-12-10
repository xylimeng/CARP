% wavecoef2logL calculates the logL contributed from a single wavelet coef
function logL = wavecoef2logL(wavecoef, level, par, sigma_hat)
nlogL1 = -log(1 - par(level + 1, 1)) + normlike([0, sigma_hat],wavecoef); % negative likelihood 
nlogL2 = -log(par(level + 1, 1)) + normlike([0, sqrt(1 + par(level + 1, 2)) * sigma_hat], wavecoef);
logL = log_expey(-nlogL1, -nlogL2);
end 
