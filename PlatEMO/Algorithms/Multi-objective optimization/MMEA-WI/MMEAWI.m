classdef MMEAWI < ALGORITHM
% <multi> <real/integer> <multimodal>
% Weighted indicator-based evolutionary algorithm for multimodal multi-objective optimization

%------------------------------- Reference --------------------------------
% W. Li, T. Zhang, R. Wang, and H. Ishibuchi, Weighted indicator-based
% evolutionary algorithm for multimodal multiobjective optimization, IEEE
% Transactions on Evolutionary Computation, 2021, 25(6): 1064-1078.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenhua Li

    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            kappa = 0.05;
            t_gen=ceil(Problem.maxFE*0.4);

            %% Generate random population
            Population = Problem.Initialization();
            [Population,pfit] = EnvironmentalSelection(Population,Problem.N,kappa,1);
            [Arcd,afit] = UpdateArc(Population,Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Arcd)
                if Problem.FE>=t_gen && rand>0.5 % Stage 2
                    [~,p1]=min(afit);
                    joint=[Arcd Population];
                    dist=pdist2(Arcd(p1).decs,joint.decs);
                    [~,so]=sort(dist,'ascend');
                    MatingPool= so(randperm(round(Problem.N/5),round(Problem.N/10))+1);
                    parents=[Arcd(p1) joint(MatingPool)];
                    Offspring  = OperatorGA(Problem,parents);
                    [Population,pfit] = EnvironmentalSelection([Population,Offspring],Problem.N,kappa,Problem.FE/Problem.maxFE);
                    [Arcd,afit] = UpdateArc(Arcd,Offspring,Problem.N);
                else % Stage 1
                    MatingPool = TournamentSelection(2,round(Problem.N/10),pfit);
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                    [Population,pfit] = EnvironmentalSelection([Population,Offspring],Problem.N,kappa,Problem.FE/Problem.maxFE);
                    [Arcd,afit] = UpdateArc(Arcd,Offspring,Problem.N);
                end
            end
        end
    end
end