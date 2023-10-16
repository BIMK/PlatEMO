function Offspring  = Operator(Problem,Population,M)
% Generate offsprings by DE and KNN based surrogate model

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = length(Population);
    
    %% Generate all candidate offsprings
    CandidateDec = OperatorDE(Problem,Population(repmat(1:N,1,M)).decs,Population(randi(N,1,N*M)).decs,Population(randi(N,1,N*M)).decs);
    
    %% Classification based preselection (CPS)
    Labels       = reshape(KNN(CandidateDec),N,M) + rand(N,M);
    [~,best]     = max(Labels,[],2);
    OffspringDec = CandidateDec((best-1)*N+(1:N)',:);
    Offspring    = Problem.Evaluation(OffspringDec);
end