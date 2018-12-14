function [MatingPool,offspringLoc] = MatingSelection(B,s)
% The mating selection of EAG-MOEA/D

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,T] = size(B);
    
    %% Select N subproblems by roulette-wheel selection
    S = sum(s,2) + 1e-6;
    D = S./sum(S) + 0.002;
    D = D./sum(D);
    offspringLoc = RouletteWheelSelection(N,1./D);
    
    %% Select two parents in each selected subproblem
    MatingPool = zeros(1,2*N);
    for i = 1 : N
        MatingPool([i,i+N]) = B(offspringLoc(i),randperm(T,2));
    end
end