% Copyright (c) 2014, J.M. Hernandez-Lobato, M.W. Hoffman, Z. Ghahramani
% This function is from the code for the paper
% Hernández-Lobato J. M., Hoffman M. W. and Ghahramani Z.
% Predictive Entropy Search for Efficient Global Optimization of Black-box
% Functions, In NIPS, 2014.
% https://bitbucket.org/jmh233/codepesnips2014
function [ ret ] = computeKnm(X, Xbar, l, sigma)

	n = size(X, 1);
	m = size(Xbar, 1);

	X = X .* repmat(sqrt(l'), n, 1);
	Xbar = Xbar .* repmat(sqrt(l'), m, 1);

	Q = repmat(sum(X.^2, 2), 1, m);
	Qbar = repmat(sum(Xbar.^2, 2)', n, 1);

	distance = Qbar + Q - 2 * X * Xbar';
	ret = sigma * exp(-0.5 * distance);
