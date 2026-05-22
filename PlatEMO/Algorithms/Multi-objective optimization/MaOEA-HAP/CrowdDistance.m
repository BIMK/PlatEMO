function MatingPool = CrowdDistance(CA,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [FrontNo,~] = NDSort(CA.objs,CA.cons,N);
    CrowdDis    = CrowdingDistance(CA.objs,FrontNo);
    for i = 1 : N
        Index = randperm(N,2);
        if CrowdDis(Index(1)) < CrowdDis(Index(2))
            MatingPool(i) = Index(2);
        elseif CrowdDis(Index(1)) > CrowdDis(Index(2))
            MatingPool(i) = Index(1);
        else
            if rand() < 0.5
                MatingPool(i) = Index(1);
            else
                MatingPool(i) = Index(2);
            end
        end
    end
end