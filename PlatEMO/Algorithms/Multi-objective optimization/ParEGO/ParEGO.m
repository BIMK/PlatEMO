classdef ParEGO < ALGORITHM
% <multi> <real/integer> <expensive>
% Efficient global optimization for Pareto optimization
% IFEs --- 10000 --- Internal GA evals per iteration

%------------------------------- Reference --------------------------------
% J. Knowles, ParEGO: A hybrid algorithm with on-line landscape
% approximation for expensive multiobjective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2006, 10(1): 50-66.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            IFEs = Algorithm.ParameterSet(10000);

            %% Generate the weight vectors and random population
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            N             = 11*Problem.D-1;
            PopDec        = UniformPoint(N,Problem.D,'Latin');
            Population    = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*PopDec+repmat(Problem.lower,N,1));
            theta         = 10.*ones(1,Problem.D);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Randomly select a weight vector and preprocess the data
                lamda  = W(randi(size(W,1)),:); 
                PopObj = Population.objs;
                [N,D]  = size(Population.decs);
                PopObj = (PopObj-repmat(min(PopObj,[],1),N,1))./repmat((max(PopObj,[],1)-min(PopObj,[],1)),N,1);
                PCheby = max(PopObj.*repmat(lamda,N,1),[],2)+0.05.*sum(PopObj.*repmat(lamda,N,1),2); 
                if N > 11*D-1+25
                    [~,index] = sort(PCheby);
                    Next      = index(1:11*D-1+25);
                else
                    Next = true(N,1);
                end
                PDec   = Population(Next).decs;
                PCheby = PCheby(Next);

                % Eliminate the solutions having duplicated inputs or outputs
                [~,distinct1] = unique(round(PDec*1e6)/1e6,'rows');
                [~,distinct2] = unique(round(PCheby*1e6)/1e6);
                distinct = intersect(distinct1,distinct2);
                PDec     = PDec(distinct,:);
                PCheby   = PCheby(distinct);

                % Surrogate-assisted prediction
                dmodel     = dacefit(PDec,PCheby,'regpoly1','corrgauss',theta,1e-5.*ones(1,D),20.*ones(1,D));
                theta      = dmodel.theta;
                PopDec     = EvolALG(Problem,PCheby,Population.decs,dmodel,IFEs);
                Population = [Population,Problem.Evaluation(PopDec)];
            end
        end
    end
end