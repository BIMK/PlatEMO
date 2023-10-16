classdef OFA < ALGORITHM
% <single> <real/integer> <large/none> <constrained/none>
% Optimal foraging algorithm

%------------------------------- Reference --------------------------------
% G. Zhu and W. Zhang, Optimal foraging algorithm for global optimization,
% Applied Soft Computing, 2017, 51: 294-313.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,rank]   = sort(FitnessSingle(Population));
            Population = Population(rank);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                PopDec     = Population.decs;
                OffDec     = PopDec + Problem.FE./Problem.maxFE.*(rand(size(PopDec))-rand(size(PopDec))).*(PopDec-PopDec([end,floor(unifrnd(ones(1,end-1),2:end))],:));
                Offspring  = Problem.Evaluation(OffDec);
                lambda     = rand(Problem.N,1);
                replace    = lambda.*FitnessSingle(Offspring)./(1+lambda*ceil(Problem.FE/Problem.N)) < FitnessSingle(Population)/ceil(Problem.FE/Problem.N);
                Population(replace) = Offspring(replace);
                [~,rank]   = sort(FitnessSingle(Population));
                Population = Population(rank);
            end
        end
    end
end