function MatingPool = MatingSelection(PopObj,K)
% The mating selection of SPEA/R

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);
    
    %% The Euclidean distance between each two solutions
    Dis = pdist2(PopObj,PopObj);
    Dis(logical(eye(N))) = inf;
    
    %% Randomly select one solution for each solution
    MatingPool = zeros(1,N);
    for i = 1 : N
        Candidates    = randperm(N,min(K,N));
        [~,nearest]   = min(Dis(i,Candidates));
        MatingPool(i) = Candidates(nearest);
    end
end