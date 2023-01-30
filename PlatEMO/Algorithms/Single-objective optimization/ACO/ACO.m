classdef ACO < ALGORITHM
% <single> <permutation> <large/none>
% Ant colony optimization

%------------------------------- Reference --------------------------------
% M. Dorigo and G. D. Caro, Ant colony optimization: a new meta-heuristic,
% Proceedings of the IEEE Congress on Evolutionary Computation, 1999.
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
            Tau = ones(Problem.D);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Move the ants
                PopDec = zeros(Problem.N,Problem.D);
                PopDec(:,1) = randi(Problem.D,Problem.N,1);
                for i = 1 : Problem.N
                    for j = 1 : Problem.D-1
                        Remain = setdiff(1:Problem.D,PopDec(i,1:j));
                        next = RouletteWheelSelection(1,Problem.C(PopDec(i,j),Remain)./Tau(PopDec(i,j),Remain));
                        PopDec(i,j+1) = Remain(next);
                    end
                end
                Population = Problem.Evaluation(PopDec);
                % Update the pheromone matrix
                dTau = zeros(Problem.D) + 1e-6;
                for i = 1 : Problem.N
                    for j = 1 : Problem.D-1
                        dTau(PopDec(i,j),PopDec(i,j+1)) = dTau(PopDec(i,j),PopDec(i,j+1)) + 1/Population(i).obj;
                    end
                    dTau(PopDec(i,end),PopDec(i,1)) = dTau(PopDec(i,end),PopDec(i,1)) + 1/Population(i).obj;
                end
                Tau = 0.5*Tau + dTau;
            end
        end
    end
end