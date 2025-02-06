function [Fitness,TDec,TMask,TempPop]  = FitnessCal(Problem,SampleNum)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    REAL    = any(Problem.encoding==1);
    TDec    = [];
    TMask   = [];
    TempPop = [];
    if REAL
        % Latin hypercube sampling
        Fitness = zeros(1,Problem.D);
        DecMat  = repmat((Problem.upper - Problem.lower),SampleNum,1).*lhsdesign(SampleNum,Problem.D) - repmat(Problem.lower,SampleNum,1);
        for i = 1 : SampleNum
            Dec  = repmat(DecMat(i,:),Problem.D,1);
            Mask = eye(Problem.D);
            Population = Problem.Evaluation(Dec.*Mask);
            TDec       = [TDec;Dec];
            TMask      = [TMask;Mask];
            TempPop    = [TempPop,Population];
            Fitness    = Fitness + sum(Population.objs,2)';
        end        
        % To reduce computation and support parallelism
        if Problem.D > 0
            AllSample = randperm(length(TempPop));
            FinalSample = AllSample(1:Problem.D);
            TempPop = TempPop(FinalSample);
            TDec = TDec(FinalSample,:);
            TMask = TMask(FinalSample,:);
        end
        Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
        Dec(:,Problem.encoding==4) = 1;
        Mask = false(Problem.N,Problem.D);
        for i = 1 : Problem.N
            Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
        end
        Population = Problem.Evaluation(Dec.*Mask);
        TDec       = [TDec;Dec];
        TMask      = [TMask;Mask];
        TempPop    = [TempPop,Population];
    else
        Fitness = zeros(1,Problem.D);
        for i = 1
            Dec = ones(Problem.D,Problem.D);
            Mask       = eye(Problem.D);
            Population = Problem.Evaluation(Dec.*Mask);
            TDec       = [TDec;Dec];
            TMask      = [TMask;Mask];
            TempPop    = [TempPop,Population];
            Fitness = Fitness + sum(Population.objs,2)';
        end
        Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
        Dec(:,Problem.encoding==4) = 1;
        Mask = false(Problem.N,Problem.D);
        for i = 1 : Problem.N
            Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
        end
        Population = Problem.Evaluation(Dec.*Mask);
        TDec       = [TDec;Dec];
        TMask      = [TMask;Mask];
        TempPop    = [TempPop,Population];
    end   
end