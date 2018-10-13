function Labels = KNN(varargin)
% K-nearest neighbor classification

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

persistent model

    if nargin == 2
        %% Train
        model.data  = [varargin{1};varargin{2}];
        model.label = [true(1,size(varargin{1},1)),false(1,size(varargin{2},1))];
    else
        %% Predict
        Distance = pdist2(varargin{1},model.data);
        [~,rank] = sort(Distance,2);
        Labels   = sum(model.label(rank(:,1:5))==1,2) > 2;
    end
end