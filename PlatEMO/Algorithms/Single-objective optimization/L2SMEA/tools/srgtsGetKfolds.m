function [IdxFolds] = srgtsGetKfolds(P, kfold, method)
%Function srgtsGetKfolds generates cross-validation sets of data following
%the "k-fold" strategy. According to this strategy, after dividing the
%available data (npoints points) into npoints/kfold clusters, each fold
%is constructed using a point randomly selected (without replacement) from
%each of the clusters. Of the k folds, a single fold is retained as the
%validation data for testing the model, and the remaining k-1 folds are
%used as training data. The cross-validation process is then repeated k
%times with each of the k folds used exactly once as the validation data.
%Note that k-fold turns out to be the leave-one out when k = npoints. It
%is important to point out that npoints/k must be an integer. Thus, for
%example:
%
%     IDXFOLDS = srgtsGetKfolds(P, K, METHOD): returns IDXFOLDS, which is the
%     vector of indexes organized such that for npoints/k clusters, as follows:
%          * Indexes for Cluster #1  : IDXFOLDS(1:k)
%          * Indexes for Cluster #2  : IDXFOLDS(k+1:2k);
%          * Indexes for Cluster #3  : IDXFOLDS(2k+1:3k);
%             ...
%          * Indexes for Cluster #n/k: IDXFOLDS((npoints-k+1:npoints);
%
%     METHOD can be:
%          * 'MaxMin', where the computation is performed as follows:
%               1. use srgtsOptimalSubDOE routine to extract npoints/kfold
%               points of the set, using 'MaxMin' criterion over 100000
%               iterations.
%               2. Remove these points from the set and go back to step 1
%               using the remaining points.
%               This does away with the clusters, and goes directly to the
%               folds. This implementation is especially successful in
%               high-dimensions, when there are no clusters anyhow.
%
%          * 'ClusterBal', which uses the balanced clustering algorithm for
%          generation of the clusters and the Anderberg's algorithm to
%          select the initial centers.
%
%Example:
%     % create a DOE
%     npoints = 200;
%     ndv = 2;
%     P = lhsdesign(npoints, ndv, ...
%                  'criterion', 'maximin', ...
%                  'iterations', 1000);
%
%     % generate the 20 k-fold sets with 10 points each
%     kfold = 20;
%     IdxFolds = srgtsGetKfolds(P, kfold, 'MaxMin');
%
%See also srgtsGetKfoldValues and srgtsComputeCrossValidation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% prepare data structures
npoints   = length( P(:,1) );
IdxFolds   = [1:kfold]';
NbClusters = npoints / kfold;

if NbClusters ~= round(NbClusters)
    msg = 'npoints/kfold must be integer.';
    error('SurrogatesToolBox:srgtsGetKfolds:',msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% run
if kfold ~= npoints

    switch method
        case 'MaxMin'
            IdxFolds = srgtsGetKfoldMaxMin(P, npoints, kfold, NbClusters);

        case 'ClusterBal'
            IdxClusters = clusterbal(P, npoints, kfold, NbClusters);

            %pick one element from each cluster
            for c1 = 1 : kfold
                idx = [c1: kfold :npoints];
                IdxFolds = vertcat(IdxFolds, IdxClusters(idx));
            end

    end
end

return
