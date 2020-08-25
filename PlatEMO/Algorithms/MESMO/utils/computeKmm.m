% Copyright (c) 2014, J.M. Hernandez-Lobato, M.W. Hoffman, Z. Ghahramani
% This function is from the code for the paper
% Hernández-Lobato J. M., Hoffman M. W. and Ghahramani Z.
% Predictive Entropy Search for Efficient Global Optimization of Black-box
% Functions, In NIPS, 2014.
% https://bitbucket.org/jmh233/codepesnips2014
function [ ret ] = computeKmm(Xbar, l, sigma, sigma0)

	m = size(Xbar, 1);
	Xbar = Xbar .* repmat(sqrt(l'), m, 1);
	Qbar = repmat(sum(Xbar.^2, 2), 1, m);
	distance = Qbar + Qbar' - 2 * (Xbar * Xbar');
	ret = sigma * exp(-0.5 * distance) + eye(m) .* (sigma0 + sigma * 1e-10);
