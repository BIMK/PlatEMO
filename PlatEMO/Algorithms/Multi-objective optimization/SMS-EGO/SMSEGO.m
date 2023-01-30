classdef SMSEGO < ALGORITHM
% <multi> <real/integer> <expensive>
% S-metric-selection-based efficient global optimization
% wmax --- 10000 --- The maximum number of internal evluation

%------------------------------- Reference --------------------------------
% W. Ponweiser, T. Wagner, D. Biermann, and M. Vincze, Multiobjective
% optimization on a limited budget of evaluations using model-assisted
% S-metric selection, Proceedings of the International Conference on
% Parallel Problem Solving from Nature, 2008, 784-794.
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
            wmax = Algorithm.ParameterSet(10000);

            %% Generate initial population based on Latin hypercube sampling
            N          = 11*Problem.D-1;
            P          = UniformPoint(N,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*P+repmat(Problem.lower,N,1));
            THETA      = 5.*ones(Problem.M,Problem.D);
            Model      = cell(1,Problem.M);
            % Lower confident bound parameter
            alpha = 1/normcdf(0.5+1/2^Problem.M);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Delete duplicated solutions
                [~,index]  = unique(Population.decs,'rows');
                Population = Population(index);
                PopDec     = Population.decs;
                PopObj     = Population.objs;

                % Optimization
                for i = 1 : Problem.M
                    [dmodel,~] = dacefit(PopDec,PopObj(1:end,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),20.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    THETA(i,:) = dmodel.theta;
                end
                w = 0;
                while w < wmax
                    drawnow('limitrate');
                    OffDec = OperatorGA(Problem,PopDec);
                    PopDec = [PopDec;OffDec]; 
                    [N,~]  = size(PopDec);
                    OffObj = zeros(N,Problem.M);
                    sigma2 = zeros(N,Problem.M);
                    for i = 1 : size(PopDec,1)
                        for j = 1 : Problem.M
                            [OffObj(i,j),~,sigma2(i,j)] = predictor(PopDec(i,:),Model{j});
                        end
                    end
                    ypot      = OffObj-alpha.*sqrt(sigma2);                           %Lower confidence bound
                    fit       = Contribution(Problem.maxFE-Problem.FE,Population(NDSort(Population.objs,1)==1).objs,ypot);
                    [~,index] = sort(fit);
                    PopDec 	  = PopDec(index(1:floor(N/2)),:);
                    fit       = fit(index(1:floor(N/2)),:);
                    w = w + floor(N/2);
                end
                [~,INDEX]  = min(fit);
                PopNew     = Problem.Evaluation(PopDec(INDEX,:));
                Population = [Population,PopNew];
            end
        end
    end
end