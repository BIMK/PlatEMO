function [MatingPool,offspringLoc] = MatingSelection(B,s)
% The mating selection of EAG-MOEA/D

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,T] = size(B);

    %% Calculate the probability of selecting each subproblem
    S = sum(s,2) + 1e-6;
    D = S./sum(S) + 0.002;
    p = D./sum(D);
    
    %% Select subproblems and parents
    % The cumulative probability
    P = cumsum(p);
    % Select N subproblems and 2 parents in each subproblem
    offspringLoc = zeros(1,N);
    MatingPool   = zeros(1,2*N);
    for i = 1 : N
        offspringLoc(i)     = find(rand<=P,1);
        MatingPool([i,i+N]) = B(offspringLoc(i),randperm(T,2));
    end
end