function Offspring  = GenerateOffsprings(Global,Population,M)
% Generate offsprings by DE and KNN based surrogate model

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = length(Population);
    
    %% Generate all candidate offsprings
    CandidateDec = Global.VariationDec(Population([repmat(1:N,1,M),randi(N,1,2*M*N)]).decs,inf,@DE);
    
    %% Classification based preselection (CPS)
    Labels       = reshape(KNN(CandidateDec),N,M) + rand(N,M);
    [~,best]     = max(Labels,[],2);
    OffspringDec = CandidateDec((best-1)*N+(1:N)',:);
    Offspring    = INDIVIDUAL(OffspringDec);
end