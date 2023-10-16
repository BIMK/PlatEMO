classdef cDPEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained dual-population evolutionary algorithm

%------------------------------- Reference --------------------------------
% M. Ming, A. Trivedi, R. Wang, and D. Srinivasan, A dual-population based
% evolutionary algorithm for constrained multi-objective optimization, IEEE
% Transactions on Evolutionary Computation, 2021, 25(4): 739-753.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Mengjun Ming

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1 = Problem.Initialization(); 
            Population2 = Population1; 
            alpha       = 2./(1+exp(1).^(-Problem.FE*10/Problem.maxFE))-1;
            para        = ceil(Problem.maxFE/Problem.N)/2 - ceil(Problem.FE/Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                % Re-rank
                Population1 = Population1(randperm(Problem.N));
                Population2 = Population2(randperm(Problem.N));
                Lia   = ismember(Population2.objs,Population1.objs, 'rows');
                gamma = 1-sum(Lia)/Problem.N;            
                [Population1,FrontNo1] = EnvironmentalSelection(Population1,Problem.N,alpha);
                [Population2,FrontNo2] = EnvironmentalSelection_noCon(Population2,Problem.N,alpha,gamma,para);

                % Offspring Reproduction
                Population_all   = [Population1,Population2];
                RankSolution_all = [FrontNo1,FrontNo2];
                MatingPool = TournamentSelection(2,2*Problem.N,RankSolution_all);
                Offspring  = OperatorGAhalf(Problem,Population_all(MatingPool));

                % Environmental Selection
                alpha = 2./(1+exp(1).^(-Problem.FE*10/Problem.maxFE)) - 1; 
                para  = ceil(Problem.maxFE/Problem.N)/2 - ceil(Problem.FE/Problem.N);
                [Population1,~] = EnvironmentalSelection([Population1,Offspring],Problem.N,alpha);
                [Population2,~] = EnvironmentalSelection_noCon([Population2,Offspring],Problem.N,alpha,gamma,para);
            end
        end
    end
end