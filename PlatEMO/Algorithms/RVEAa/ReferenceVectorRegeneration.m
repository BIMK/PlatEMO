function V = ReferenceVectorRegeneration(PopObj,V)
% Reference vector regeneration

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj        = PopObj - repmat(min(PopObj,[],1),size(PopObj,1),1);
    [~,associate] = max(1-pdist2(PopObj,V,'cosine'),[],2);
    inValid       = setdiff(1:size(V,1),associate);
    V(inValid,:)  = rand(length(inValid),size(V,2)).*repmat(max(PopObj,[],1),length(inValid),1);
end