function [Pbest,Gbest] = GetBest(Archive,W,Z)
% Update the pbest and gbest of each particle

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select the pbest
    AObj  = Archive.objs - repmat(Z,length(Archive),1);
    NormW = sqrt(sum(W.^2,2))';
    d1    = AObj*W'./repmat(NormW,length(Archive),1);
    d2    = sqrt(repmat(sum(AObj.^2,2),1,size(W,1))-d1.^2);
    PBI   = d1 + 5*d2;
    [~,p] = min(PBI,[],1);
    Pbest = Archive(p);
    
    %% Select the gbest
    Gbest = Archive(randi(size(Archive,1),1,size(W,1)));
end