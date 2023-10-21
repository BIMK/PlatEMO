classdef RGA_M1_2 < ALGORITHM
% <multi> <real> <expensive> <constrained>
% Real-coded genetic algorithm with framework M1-2
% wmax --- 20 --- Number of generations before updating surrogate models
% mu   ---  5  -- Number of real evaluated solutions at each iteration

%------------------------------- Reference --------------------------------
% K. Deb, R. Hussein, P. C. Roy, and G. Toscano-Pulido, A taxonomy for 
% metamodeling framework for evolutionary multiobjective optimization, IEEE
% Transactons on Evolutionary Computation, 2019, 23(1): 104-116.
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
            %% Parameter settings
            [wmax,mu]   = Algorithm.ParameterSet(20,5);
            
            %% Generate the initial population
            NI          = 11*Problem.D-1;
            P           = lhsamp(NI,Problem.D);
            Population  = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            Archive     = UpdateArchive(Population,Problem.N);
            THETA       = 5.*ones(Problem.M+size(Population.cons,2),Problem.D);
            
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                PopDec  = Population.decs;
                PopCon  = max(0,Population.cons);
                PopCon  = (PopCon-min(PopCon,[],1))./(max(PopCon,[],1)-min(PopCon,[],1)); 
                PopCon(:,isnan(PopCon(1,:))) = 0;
                PopObj  = [Population.objs,PopCon];
                M       = size(PopObj,2);
                Model   = cell(1,M);
                
                %% Construct m+q surrogate models for each objectives and constraints
                for i = 1 : M
                    dmodel     = dacefit(PopDec,PopObj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    THETA(i,:) = dmodel.theta;
                end
                
                %% Use the surrogates
                [FrontNo,~] = NDSort(Population.objs,sum(PopCon,2),inf);
                CrowdDis    = CrowdingDistance(Population.objs,FrontNo);
                for w = 1:wmax
                    drawnow();
                    MatingPool = TournamentSelection(2,NI,FrontNo,-CrowdDis);
                    OffDec     = OperatorGA(Problem,PopDec(MatingPool,:));
                    PopDec     = cat(1,PopDec,OffDec);
                    [N,~]      = size(PopDec);
                    PopObj     = zeros(N,M);
                    MSE        = zeros(N,M);
                    for i = 1: N
                        for j = 1 : M
                            [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),Model{j});
                        end
                    end
                    [PopDec,PopObj,FrontNo,CrowdDis] = EnvironmentalSelection(PopDec,PopObj,Problem.M,NI);
                end
                
                %% Select mu solutions for real evaluation
                [NewDec,~,~,~] = EnvironmentalSelection(PopDec,PopObj,Problem.M,mu);
                New            = Problem.Evaluation(NewDec);
                
                %% Update population P and archive A
                Population  = UpdatePopulation(Population,New,NI-mu);
                Archive     = cat(2,Archive,New);
                Archive     = UpdateArchive(Archive,Problem.N);
            end
        end
    end
end