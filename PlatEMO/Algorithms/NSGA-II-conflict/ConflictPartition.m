function Psi = ConflictPartition(PopObj,NS)
% Partitioning of objectives using conflict information

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    M = size(PopObj,2);
    
    %% Calculate the modified correlation matrix of objectives
    cMatrix = 1 - corr(PopObj);
    cMatrix = repmat(max(cMatrix,[],2)',M,1).*(1-eye(M));
    
    %% Partition the objectives
    k      = ceil(M/NS);
    Psi    = cell(1,NS);
    remain = 1 : M;
    for s = 1 : NS-1
        [~,L]  = sort(cMatrix(remain,remain),2);
        [~,i]  = min(L(:,k));
        Psi{s} = remain(L(i,1:k));
        remain(L(i,1:k)) = [];
    end
    Psi{NS} = remain;
end