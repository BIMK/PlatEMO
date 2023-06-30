classdef ADSAPSO < ALGORITHM
% <multi/many> <real/integer> <expensive>
% Adaptive dropout based surrogate-assisted particle swarm optimization
% k    ---   5 --- Number of re-evaluated solutions
% beta --- 0.5 --- Percentage of Dropout

%------------------------------- Reference --------------------------------
% J. Lin, C. He, and R. Cheng, Adaptive dropout for high-dimensional 
% expensive multiobjective optimization, Complex & Intelligent Systems,
% 2022, 8(1): 271¨C285.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jianqing Lin

    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            [k,beta] = Algorithm.ParameterSet(5,0.5);
            Init_Num = 100;     % Initial number of solutions
            N_a      = 200;     % The number of solutions for building surrogate models
            N_s      = 50;      % The number of the well- and poorly performing solutions

            %% Generate initial population
            InitDec  = repmat((Problem.upper - Problem.lower),Init_Num, 1).* lhsdesign(Init_Num, Problem.D) + repmat(Problem.lower, Init_Num, 1);  % lhs design
            Arc      = Problem.Evaluation(InitDec);
            
            %% Optimization
            while Algorithm.NotTerminated(Arc)
                Offspring  = Operator(Problem,Arc,k,beta,N_a,N_s);
                Offspring  = Problem.Evaluation(Offspring);
                Arc        = [Arc,Offspring];
            end
        end
    end
end