function Density = DensityEstimate(Population,W)
% Estimate the density of each individual

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    Zmax       = max(Population.objs,[],1);
    Zmin       = min(Population.objs,[],1);
    SPopObj    = (Population.objs-repmat(Zmin,size(Population.objs,1),1))./(repmat(Zmax,size(Population.objs,1),1)-repmat(Zmin,size(Population.objs,1),1));
    [~,Region] = max(1-pdist2(SPopObj,W,'cosine'),[],2);
    [value,~]  = sort(Region,'ascend');
    flag       = max(value);
    counter    = histc(value,1:flag);
    Density    = counter(Region);
end
