classdef DirHVEI < ALGORITHM
% <2024> <multi/many> <real/integer> <expensive>
% Expected direction-based hypervolume improvement
% batch_size --- 5 --- number of true function evaluations per iteration

%------------------------------- Reference --------------------------------
% L. Zhao and Q. Zhang. Hypervolume-guided decomposition for parallel
% expensive multiobjective optimization. IEEE Transactions on Evolutionary
% Computation, 28(2): 432-444, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function was written by Liang Zhao.
% https://github.com/mobo-d/DirHV-EGO

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            batch_size = Algorithm.ParameterSet(5); 
            % number of initial samples
            n_init = 11*Problem.D-1;
            % Initial hyperparameters for GP
            theta = repmat({(n_init ^ (-1 ./ n_init)) .* ones(1, Problem.D)}, 1, Problem.M);
            
            %% Generate initial design using LHS or other DOE methods
            x_lhs   = lhsdesign(n_init, Problem.D,'criterion','maximin','iterations',1000);
            x_init  = Problem.lower +  (Problem.upper - Problem.lower).*x_lhs;  
            Archive = Problem.Evaluation(x_init);     
            % find non-dominated solutions
            FrontNo = NDSort(Archive.objs,1); 
			
            %% Optimization
            while Algorithm.NotTerminated(Archive(FrontNo==1))
              %% Scale the objective values 
                train_x = Archive.decs; ori_objs = Archive.objs;
                ymin    = min(ori_objs,[],1); ymax = max(ori_objs,[],1);
		        train_y = (ori_objs-ymin)./(ymax - ymin);
                train_y_nds = train_y(FrontNo==1,:);

              %% Bulid GP model for each objective function 
                GPModels = cell(1,Problem.M);
                for i = 1 : Problem.M
                    GPModels{i} = Dacefit(train_x,train_y(:,i),'regpoly0','corrgauss',theta{i},1e-6*ones(1,Problem.D),20*ones(1,Problem.D));
                    theta{i}    = GPModels{i}.theta;
                end 
                
              %% Maximize DirHV-EI using the MOEA/D framework and select multiple candidate points
                Batch_size = min(Problem.maxFE - Problem.FE,batch_size); % the total budget is  Problem.maxFE 
                new_x = Opt_DirHV_EI(Problem.M,Problem.D,Problem.lower,Problem.upper,GPModels,train_y_nds,Batch_size);  
              
              %% Expensive Evaluation
                Archive = [Archive,Problem.Evaluation(new_x)];
                FrontNo = NDSort(Archive.objs,1);
            end
        end
    end
end