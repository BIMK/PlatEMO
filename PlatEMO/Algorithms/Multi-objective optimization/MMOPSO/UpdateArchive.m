function A = UpdateArchive(A,N)
% Update the external archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select only the non-dominated solutions in the archive
    A = A(NDSort(A.objs,1)==1);
    
    %% Sort the solutions in A according to their crowding distances
    [~,rank] = sort(CrowdingDistance(A.objs),'descend');
    A        = A(rank(1:min(N,end)));
end