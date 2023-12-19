function [FitnessInit,FitnessOpt,SparseRate,TDec,TMask,TempPop] = FitnessCal(Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
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
        Fitness  = zeros(5,Problem.D);
        Interval = (Problem.upper-Problem.lower)./5;
        for i = 1 : 5
            for j = 1 : 2
                Dec        = unifrnd(repmat(Problem.lower + Interval*(i-1),Problem.D,1),repmat(Problem.lower + Interval*(i),Problem.D,1));
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                Fitness(i,:) =    Fitness(i,:) + NDSort([Population.objs,Population.cons],inf);
            end        
        end
        % To reduce computation and support parallelism
        if Problem.D > 2000
            AllSample   = randperm(length(TempPop));
            FinalSample = AllSample(1:Problem.D);
            TempPop     = TempPop(FinalSample);
            TDec        = TDec(FinalSample,:);
            TMask       = TMask(FinalSample,:);
        end
        FitnessInit = sum(Fitness);
        FitnessOpt  = kmeans(sum(Fitness)',2); 
        [SparseRate,FitnessOpt] = ClusterLabel(FitnessInit,FitnessOpt);
        FitnessOpt = FitnessOpt';
    else
        Fitness = zeros(1,Problem.D);
        for i = 1
            Dec          = ones(Problem.D,Problem.D);
            Mask         = eye(Problem.D);
            Population   = Problem.Evaluation(Dec.*Mask);
        	TDec         = [TDec;Dec];
          	TMask        = [TMask;Mask];
           	TempPop      = [TempPop,Population];
           	Fitness(i,:) = Fitness(i,:) + NDSort([Population.objs,Population.cons],inf);
        end        
        FitnessInit = Fitness;
        FitnessOpt  = kmeans((Fitness)',2); 
        [SparseRate,FitnessOpt] = ClusterLabel(FitnessInit,FitnessOpt);
        FitnessOpt  = FitnessOpt';
        Mask        = zeros(Problem.N,Problem.D);
        Dec         = ones(Problem.N,Problem.D);
        for i = 1 : Problem.N
            Mask(i,TournamentSelection(2,ceil(rand*Problem.D),FitnessInit)) = 1;
        end           
        Population = Problem.Evaluation(Dec.*Mask);
        TDec       = [TDec;Dec];
        TMask      = [TMask;Mask];
        TempPop    = [TempPop,Population];
    end   
end

function [SparseRate,Fitness] = ClusterLabel(FitnessInit,Fitness)
    Num1 = sum(Fitness == 1);
    Num2 = sum(Fitness == 2);
    Num1Value = sum(FitnessInit(Fitness == 1));
    Num2Value = sum(FitnessInit(Fitness == 2));
    if Num1Value < Num2Value
        SparseRate = Num1/(Num1 + Num2);
        Fitness(Fitness == 1) = 11;
        Fitness(Fitness == 2) = 12;
    else   
        SparseRate = Num2/(Num1 + Num2);
        Fitness(Fitness == 1) = 12;
        Fitness(Fitness == 2) = 11;
    end
end