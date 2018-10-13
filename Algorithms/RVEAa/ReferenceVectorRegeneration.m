function V = ReferenceVectorRegeneration(PopObj,V)
% Reference vector regeneration

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj        = PopObj - repmat(min(PopObj,[],1),size(PopObj,1),1);
    [~,associate] = max(1-pdist2(PopObj,V,'cosine'),[],2);
    inValid       = setdiff(1:size(V,1),associate);
    V(inValid,:)  = rand(length(inValid),size(V,2)).*repmat(max(PopObj,[],1),length(inValid),1);
end