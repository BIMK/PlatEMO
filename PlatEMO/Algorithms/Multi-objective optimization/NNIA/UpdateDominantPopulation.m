function D = UpdateDominantPopulation(D,N)
% Update the dominant population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    D = D(NDSort(D.objs,1)==1);
    [~,rank] = sort(CrowdingDistance(D.objs),'descend');
    D = D(rank(1:min(N,length(rank))));
end