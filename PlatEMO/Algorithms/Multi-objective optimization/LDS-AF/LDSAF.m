classdef LDSAF < ALGORITHM
% <2024> <multi> <real/integer> <large/none> <expensive>
% Low-dimensional surrogate aggregation function

%------------------------------- Reference --------------------------------
% H. Gu, H. Wang, C. He, B. Yuan, and Y. Jin. Large-scale multiobjective
% evolutionary algorithm guided by low-dimensional surrogates of
% scalarization functions. Evolutionary Computation, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Haoran Gu

    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            [delta,N_s] = Algorithm.ParameterSet(0.9, 20);
            
            %% Generate the weight vectors
            [W, Problem.N] = UniformPoint(Problem.N, Problem.M);
            
            %% Detect the neighbours of each solution
            T      = ceil(Problem.N / 10);
            B      = pdist2(W, W);
            [~, B] = sort(B, 2);
            B      = B(:, 1 : T);
            
            %% Initialize population
            PopDec     = UniformPoint(Problem.N, Problem.D, 'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper - Problem.lower, Problem.N, 1) .* PopDec + repmat(Problem.lower, Problem.N, 1));
            A          = Population;
            Z          = min(Population.objs, [], 1);
            
            %% Define RBFN
            RBFN_list1 = RBFN_PP(Problem);
            RBFN_list2 = RBFN_P(Problem);
            P_center   = B;
            
            %% Optimization
            while Algorithm.NotTerminated(A)
                % For each sub-problem
                for i = 1 : Problem.N
                    %% Model-construction
                    if Problem.FE >= 350
                        pdb  = pdist2(W, W);
                        pds2 = pdb(i,B(i,:));
                        farv = find(pds2>=(max(pds2)-0.00001));
                        a    = randperm(length(farv));
                        
                        center = P_center(i,farv(a(1)));
                        
                        RBFN_list2(i) = RBFN_list2(i).ModelConstruction(A, B(i, :), W, Z,center);
                        model = RBFN_list2(i);
                    else
                        RBFN_list1(i) = RBFN_list1(i).ModelConstruction(A, B(i, :), W, Z);
                        model = RBFN_list1(i);
                    end
                    
                    %% Choose the parents
                    if rand < delta
                        P = B(i, randperm(end));
                    else
                        P = randperm(Problem.N);
                    end
                    
                    %% Solution-generation
                    y_i = SolutionSelection(Problem, Population, P, model, N_s, i);
                    
                    %% Evaluate offspring
                    y_i = Problem.Evaluation(y_i);
                    
                    %% Update the reference point
                    Z = min(Z, y_i.obj);
                    if Problem.FE <= 2*Problem.N
                        nr = 2;
                    else
                        nr = 4;
                    end
                    
                    %% Update population and archive
                    g_old = max(abs(Population(P).objs - repmat(Z, length(P), 1)) .* W(P, :), [], 2);
                    g_new = max(repmat(abs(y_i.obj - Z), length(P), 1) .* W(P, :), [], 2);
                    Population(P(find(g_old >= g_new, nr))) = y_i;
                    A = [A, y_i];
                    
                    %% Check termination criteria
                    Algorithm.NotTerminated(A);
                end
            end
        end
    end
end