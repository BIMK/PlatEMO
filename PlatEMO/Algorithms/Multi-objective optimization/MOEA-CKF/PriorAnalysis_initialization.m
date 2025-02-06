function [TDec,TMask,TempPop,Fitness]=PriorAnalysis_initialization(Problem,REAL)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    TDec    = [];
    TMask   = [];
    TempPop = [];
    Fitness = zeros(1,Problem.D);
    for i = 1 : 1+4*REAL
        if REAL
            Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
        else
            Dec = ones(Problem.D,Problem.D);
        end
        Mask       = eye(Problem.D);
        Population = Problem.Evaluation(Dec.*Mask);
        TDec       = [TDec;Dec];
        TMask      = [TMask;Mask];
        TempPop    = [TempPop,Population];
        Fitness    = Fitness + NDSort([Population.objs,Population.cons],inf);
    end 
    
    % Generate initial population
    if REAL
        Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
    else
        Dec = ones(Problem.N,Problem.D);
    end
    Mask = zeros(Problem.N,Problem.D);
    for i = 1 : Problem.N
        Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
    end
    Population = Problem.Evaluation(Dec.*Mask);
    TDec       = [TDec;Dec];
    TMask      = [TMask;Mask];
    TempPop    = [TempPop,Population];
end