function [CrowdDis]=Crowding(Pop)
%Harmonic average distance of each solution in the decision space
%Return: the crowding distance of each individual

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenhua Li

    [N, ~]=size(Pop);
    K=N-1;
    Z = min(Pop,[],1);
    Zmax = max(Pop,[],1);
    pop=(Pop-repmat(Z,N,1))./repmat(Zmax-Z,N,1);
    distance=pdist2(pop,pop);
    [value,~]=sort(distance,2);
    CrowdDis=K./sum(1./value(:,2:N),2);
end