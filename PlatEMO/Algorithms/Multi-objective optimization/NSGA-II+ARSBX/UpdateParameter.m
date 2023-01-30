function [B,m,ps,Population] = UpdateParameter(Problem,Population)
% Update the parameters in ARSBX

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopDec = Population.decs;
    Flag   = Population.adds(zeros(length(Population),1));
    Ori    = sum(Flag == 1);
    Eig    = sum(Flag == 2);
    ps     = 1/(1+exp(-Problem.M*sqrt(Problem.D)*((Ori+1)/(Eig+Ori+2)-0.5)*Problem.FE/Problem.maxFE));
    C        = cov(PopDec);
    [B,E]    = eig(C);
    E        = diag(E);
    E        = sqrt(E);
    [~,Rank] = sort(E,'descend');
    B        = B(:,Rank);
    m        = mean(PopDec);
end