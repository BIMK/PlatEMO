function [num,active] = NoActive(PopObj,V)
% Detect inactive reference vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    [N,~] = size(PopObj);
    NV    = size(V,1);
    
    %% Translate the population
    PopObj = PopObj - repmat(min(PopObj,[],1),N,1);

    %% Associate each solution to a reference vector
    Angle   = acos(1-pdist2(PopObj,V,'cosine'));
    [~,associate] = min(Angle,[],2);
    active  = unique(associate);
	num     = NV-length(active);
end