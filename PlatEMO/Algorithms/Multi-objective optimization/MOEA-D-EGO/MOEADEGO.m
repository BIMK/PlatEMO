classdef MOEADEGO < ALGORITHM
% <2010> <multi> <real/integer> <expensive>
% MOEA/D with efficient global optimization
% batch_size --- 5 --- Number of true function evaluations per iteration

%------------------------------- Reference --------------------------------
% Q. Zhang, W. Liu, E. Tsang, and B. Virginas. Expensive multiobjective
% optimization by MOEA/D with Gaussian process model. IEEE Transactions on
% Evolutionary Computation, 2010, 14(3): 456-474.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function was written by Liang Zhao (liazhao5-c@my.cityu.edu.hk).
% - The Java Code of MOEA/D-EGO (written by Wudong Liu) is avaliable at 
%   https://sites.google.com/view/moead/resources 
% - The Matlab Code of MOEA/D-EGO without FuzzyCM (written by Liang ZHAO) is 
%   avaliable at https://github.com/mobo-d/MOEAD-EGO

    methods
        function main(Algorithm,Problem)
           %% Parameter setting
            batch_size = Algorithm.ParameterSet(5); 
            
            %% Generate initial design using LHS or other DOE methods
            x_lhs   = lhsdesign(11*Problem.D-1, Problem.D,'criterion','maximin','iterations',1000);
            x_init  = Problem.lower +  (Problem.upper - Problem.lower).*x_lhs;  
            Archive = Problem.Evaluation(x_init);     
            % find non-dominated solutions
            FrontNo = NDSort(Archive.objs,1); 
           
            %% Optimization
            while Algorithm.NotTerminated(Archive(FrontNo==1))     
              %% Maximize ETI using MOEA/D and select q candidate points
                Batch_size = min(Problem.maxFE - Problem.FE,batch_size); % the total budget is  Problem.maxFE 
                train_x    = Archive.decs; train_y = Archive.objs;
                new_x      = Opt_ETI_FCM(Problem.M,Problem.D,Problem.lower,Problem.upper,Batch_size,train_x,train_y);              
 
               %% Expensive Evaluation
                Archive = [Archive,Problem.Evaluation(new_x)];
                FrontNo = NDSort(Archive.objs,1);
            end
        end
    end
end