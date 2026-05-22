classdef MOEAIB < ALGORITHM
% <2026> <multi/many> <real> <large>
% MOEA with information bottleneck
% nSel      --- 1 --- Number of selected solutions for decision variable clustering
% nPer      --- 4 --- Number of perturbations on each solution for decision variable clustering
% optimizer --- 1 --- Option of basic optimizer, WOF: 1, ReMO: 2

%------------------------------- Reference --------------------------------
% L. Si, Z. Wang, X. Zhang, and Y. Tian. Information bottleneck theory-
% guided dimension reduction for large-scale multi-objective optimization.
% IEEE Transactions on Evolutionary Computation, 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [nSel,nPer,optimizer] = Algorithm.ParameterSet(1,4,1);

            %% Generate random population
            ObjBool    = -ones(Problem.maxFE, Problem.M);
            Population = Problem.Initialization();
            ObjBool(1:Problem.N,:) = Population.objs;
            I = 0;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Problem.FE < 0.6*Problem.maxFE
                    %% Detect the importance level of each variable
                    if mod(I,10) == 0
                        [AlphaCon,AlphaDiv,ConVars,DivVars] = ImportanceLevel(Problem,Population,nSel,nPer);
                        ConDim = max(ceil(numel(ConVars)*AlphaCon),1);
                        DivDim = max(ceil(numel(DivVars)*AlphaDiv),1);
                    end
                    if optimizer == 1
                        % Optimization in the original space
                        [Population,ObjBool] = WOF_NSGAIII(Problem,Population,5,ObjBool);

                        % Optimization in the reduced space
                        [Population,ObjBool] = WOF_Optimizer(Problem,Population,ConDim,DivDim,ConVars,DivVars,ObjBool);
                    elseif optimizer == 2
                        % Optimization in the reduced space
                        [RePopulation,LowDecs,MapCon,MapDiv,LB,UB] = ReMO_AdjustMap(Problem,Population.decs,ConVars,DivVars,ConDim,DivDim);
                        Population = ReMO_Optimizer(Problem,RePopulation,Population,LowDecs,MapCon,MapDiv,ConVars,DivVars,LB,UB);
                    else
                        % Optimization in the reduced space
                        Population = MOEAPSL_Optimizer(Problem,Population,ConDim,DivDim,ConVars,DivVars);
                    end
                    I = I + 1;
                else
                    Iter = ceil((Problem.maxFE-Problem.FE)/Problem.N);
                    if optimizer == 1
                        [Population,ObjBool] = WOF_NSGAIII(Problem,Population,Iter,ObjBool);
                    else
                        Population = ReMO_NSGAII(Problem,Population,Iter);
                    end
                end
            end
        end
    end
end