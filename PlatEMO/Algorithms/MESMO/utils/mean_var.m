% Copyright (c) 2017 Zi Wang
function [meanVector, varVector, meangrad, vargrad] = mean_var(x, ...
    Xsamples, Ysamples, KernelMatrixInv, l_, sigma_, sigma0_)
% This function returns the posterior GP mean and variance predictions made
% at queried point x conditioned on observations Xsamples (size T x d), 
% Ysamples (size T x 1).
% If nargout <=2, x is of size n x d, meanVector is of size n x 1,
% and varVector is of size n x 1. 
% If nargout > 2, x must be of size 1 x d, and all outputs are scalers.
n = size(x,1);
kstar = computeKnm(x, Xsamples, l_', sigma_);
meanVector = kstar * KernelMatrixInv * Ysamples;
varVector = ones(n,1)*(sigma_ + sigma0_) - sum((kstar * KernelMatrixInv) .* kstar, 2);

if nargout > 2
    dkstar = -repmat(l_, size(Xsamples, 1), 1) .* ...
        (repmat(x, size(Xsamples, 1), 1) - Xsamples) .* ...
        repmat(kstar', 1, size(Xsamples, 2));
    meangrad = dkstar' * KernelMatrixInv * Ysamples;
    vargrad = -2*sum( (dkstar' * KernelMatrixInv) *kstar', 2);
end