classdef EIMEGO < ALGORITHM
% <multi> <real/integer> <expensive>
% Expected improvement matrix based efficient global optimization
% InfillCriterionIndex --- 1 --- infill criterion index number

%------------------------------- Reference --------------------------------
% D. Zhan, Y. Cheng, and J. Liu, Expected improvement matrix-based infill
% criteria for expensive multiobjective optimization, IEEE Transactions on
% Evolutionary Computation, 2017, 21(6): 956-975.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Dawei Zhan

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            % 1 for Euclidean distance-based EIM criterion (default)
            % 2 for Maximin distance-based EIM criterion
            % 3 for Hypervolume-based EIM criterion
            InfillCriterionIndex = Algorithm.ParameterSet(1);

            %% Generate the initial design points
            % number of design variables
            D = Problem.D;
            % number of objective functions
            M = Problem.M;
            % number of initial design points
            N  = 11*D-1;
            % generate initial design points using Latin Hypercube sampling
            PopDec = repmat(Problem.upper-Problem.lower,N,1).*UniformPoint(N,D,'Latin')+repmat(Problem.lower,N,1);
            % calculate initial design points
            Population   = Problem.Evaluation(PopDec);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % scale the objective values to 0 and 1
                PopDec = Population.decs;
                PopObj = Population.objs;
                N  = size(PopDec,1);       
                PopObjScaled = (PopObj-repmat(min(PopObj),N,1))./(repmat(max(PopObj),N,1)-repmat(min(PopObj),N,1));   
                % bulid Kriging models for all the objective functions
                KrigingModel = cell(1,M);
                for i = 1 : M
                    KrigingModel{i}= dacefit(PopDec,PopObjScaled(:,i),'regpoly0','corrgauss',1*ones(1,D),0.001*ones(1,D),1000*ones(1,D));
                end       
                % select one candidate with the maximum EIM value using GA
                PopDec     = InfillSamplingEIM(Problem,KrigingModel,PopObjScaled,InfillCriterionIndex);
                Population = [Population,Problem.Evaluation(PopDec)];
            end
        end
    end
end