classdef MGSAEA < ALGORITHM
% <multi> <real> <expensive> <constrained>
% Multigranularity surrogate-assisted constrained evolutionary algorithm
% wmax   ---    20 --- Number of generations before updating surrogate models
% mu     ---     5 --- Number of real evaluated solutions at each iteration
% gap    ---    20 --- Parameter calculating the change rate of ideal points
% lambda --- 0.001 --- Parameter determining the evolving stages 

%------------------------------- Reference --------------------------------
% Y. Zhang, H. Jiang, Y. Tian, H. Ma, and X. Zhang, Multigranularity
% surrogate modeling for evolutionary multiobjective optimization with
% expensive constraints, IEEE Transactions on Neural Networks and Learning
% Systems, 2023.
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
             %% Parameter setting
            [wmax,mu,gap,lambda] = Algorithm.ParameterSet(20,5,20,1e-3);
            
            %% Generate the initial population 
            NI          = 11*Problem.D-1;
            P           = lhsamp(NI,Problem.D);
            Population  = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));                      
            flag        = 0;  % 0 for the first stage, 1 for the second stage
            Iter        = 1;
            IdealPoints = [];
            Archive     = UpdateArchive(Population,Problem.N);
            THETA_OBJ   = 5.*ones(Problem.M,Problem.D);
            THETA_CV    = 5.*ones(1,Problem.D);
            THETA_CON   = 5.*ones(size(Population.cons,2),Problem.D);
            
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                IdealPoints(Iter,:) = min(Population.objs,[],1);
                if Iter > gap && flag == 0
                    max_change = calc_maxchange(IdealPoints,Iter,gap);
                    if max_change <= lambda
                        flag = 1;
                    end
                end
                if flag ==0 
                    PopDec  = Population.decs;
                    PopObj  = Population.objs;
                    M       = Problem.M;
                    Model   = cell(1,M);
                    % Construct M surrogate models for M objectives
                    for i = 1 : M
                        dmodel     = dacefit(PopDec,PopObj(:,i),'regpoly0','corrgauss',THETA_OBJ(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        Model{i}   = dmodel;
                        THETA_OBJ(i,:) = dmodel.theta;
                    end
                    % Use the surrogates
                    Fitness = CalFitness(PopObj);
                    for w = 1:wmax
                        drawnow();
                        MatingPool = TournamentSelection(2,NI,Fitness);
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
                        [PopDec,PopObj,Fitness] = EnvironmentalSelection(PopDec,PopObj,NI,Problem.M);
                    end
                    % Select mu solutions for real evaluation in the first stage
                    [NewDec,~,~] = EnvironmentalSelection(PopDec,PopObj,mu,Problem.M);
                    New          = Problem.Evaluation(NewDec);
                    % Update Population and Archive with mu new solutions
                    Population   = UpdatePopulation(Population,New,NI-mu);
                    Archive      = cat(2,Archive,New);
                    Archive      = UpdateArchive(Archive,Problem.N);
                else  
                    MaxCon = max(max(0,Population.cons));
                    Ninf   = length(find(MaxCon>0)); % The number of not fully staisfied constraints
                    if Ninf == size(Population.cons,2)
                        status  = 1;
                        CV      = NormalizeCV(Population.cons); 
                        PopDec  = Population.decs;
                        PopObj  = [Population.objs,CV];
                        M       = Problem.M+1;
                        Model   = cell(1,M);
                        Fitness = CalFitness(Population.objs,Population.cons);
                        THETA   = [THETA_OBJ;THETA_CV];
                        % Construct M+1 surrogate models for M objectives and the CV function
                        for i = 1 : M
                            dmodel     = dacefit(PopDec,PopObj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                            Model{i}   = dmodel;
                            THETA(i,:) = dmodel.theta;
                        end
                        THETA_OBJ = THETA(1:Problem.M,:);
                        THETA_CV  = THETA(end,:);
                    elseif Ninf > 0 && Ninf < size(Population.cons,2)
                        status   = 2;
                        PopCon   = max(0,Population.cons);
                        MaxCon   = max(PopCon);
                        index    = find(MaxCon>0);
                        PopCon(:,MaxCon==0) = [];
                        PopCon   = (PopCon-min(PopCon,[],1))./(max(PopCon,[],1)-min(PopCon,[],1));
                        PopDec   = Population.decs;
                        PopObj   = [Population.objs,PopCon];
                        M        = size(PopObj,2);
                        Model    = cell(1,M);
                        Fitness  = CalFitness([Population.objs,sum(max(0,Population.cons),2)]);
                        THETA    = [THETA_OBJ;THETA_CON(index,:)];
                        % Construct surrogate models for M objectives and the not
                        % fully satisfied constraints
                        for i = 1 : M
                            dmodel     = dacefit(PopDec,PopObj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                            Model{i}   = dmodel;
                            THETA(i,:) = dmodel.theta;
                        end
                        THETA_OBJ = THETA(1:Problem.M,:);
                        THETA_CON(index,:) = THETA(Problem.M+1:end,:);
                    elseif Ninf  == 0
                        status   = 3;
                        PopDec   = Population.decs;
                        PopObj   = Population.objs;
                        M        = Problem.M;
                        Model    = cell(1,M);
                        Fitness  = CalFitness(PopObj);
                        % Construct M surrogate models for M objectives
                        for i = 1 : M
                            dmodel     = dacefit(PopDec,PopObj(:,i),'regpoly0','corrgauss',THETA_OBJ(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                            Model{i}   = dmodel;
                            THETA_OBJ(i,:) = dmodel.theta;
                        end
                    end
                    % Use the surrogates
                    for w = 1: wmax
                        drawnow();
                        MatingPool = TournamentSelection(2,NI,Fitness);
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
                        [PopDec,PopObj,Fitness] = EnvironmentalSelection(PopDec,PopObj,NI,Problem.M,status);
                    end
                    % Select mu solutions for real evaluation in the second stage
                    [NewDec,~,~] = EnvironmentalSelection(PopDec,PopObj,mu,Problem.M,status);
                    New          = Problem.Evaluation(NewDec);
                    % Update Population and Archive with mu new solutions
                    Population   = UpdatePopulation(Population,New,NI-mu,status);
                    Archive      = cat(2,Archive,New);
                    Archive      = UpdateArchive(Archive,Problem.N);
                end
                Iter = Iter +1;
            end
        end
    end
end

function max_change = calc_maxchange(ideal_points,Iter,gap)
    % Calculate the maximum change rate of ideal points
    delta = 1e-6 * ones(1,size(ideal_points,2));
    rz    = abs((ideal_points(Iter,:) - ideal_points(Iter - gap,:)) ./ max(ideal_points(Iter - gap,:),delta));  
    max_change = max(rz);
end

function CV = NormalizeCV(PopCon)
    % Calculate the normalized overall constraints violation
    PopCon = max(0,PopCon);
    PopCon = (PopCon-min(PopCon,[],1))./(max(PopCon,[],1)-min(PopCon,[],1));
    PopCon(:,isnan(PopCon(1,:))) = 0;
    CV = sum(PopCon,2);
end