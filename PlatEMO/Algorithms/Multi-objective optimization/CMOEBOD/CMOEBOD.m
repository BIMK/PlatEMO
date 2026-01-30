classdef CMOEBOD < ALGORITHM
% <2025> <multi/many> <real> <expensive> <constrained>
% Constrained multiobjective evolutionary Bayesian optimization based on decomposition
% wmax --- 20 --- Number of generations before updating Kriging models
% mu   ---  5 --- Number of re-evaluated solutions at each generation

%------------------------------- Reference --------------------------------
% Z. Zhang, Y. Wang, G.Sun, and T. Pang. A novel evolutionary Bayesian 
% optimization algorithm based on decomposition for expensive constrained 
% multiobjective optimization problems. IEEE Transactions on Systems, Man, 
% and Cybernetics: Systems, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao.zhang.cn@gmail.com)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [wmax,mu] = Algorithm.ParameterSet(20,5);

            %% Generate the reference points and population
            global V beta NI 
            [V,Problem.N] = UniformPoint(Problem.N,Problem.M);
            NI            = Problem.N;
            P             = UniformPoint(NI,Problem.D,'Latin');
            DB            = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            Population    = DB;           
            THETA_obj     = 5.*ones(Problem.M,Problem.D);
            THETA_con     = 5.*ones(size(DB.cons,2),Problem.D);
            beta          = 0.1;

            %% Optimization
            while Algorithm.NotTerminated(DB)
                % Model Construction
                [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(DB,THETA_obj,THETA_con);
                
                % Collaborative Evolutionary Optimization
                P = optimizaiton(Population,wmax,Model_obj,Model_con,Problem);
                
                % Bilevel Candidate Selection
                DB = NewSelect(P,DB,mu,Problem);   
                
                % Update Population
                Population = PUpdate(DB,NI);           
            end
        end
    end
end

function [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(A2,THETA_obj,THETA_con)
    Dec = A2.decs;
    Obj = A2.objs;
    Con = A2.cons;
    Len_dec = size(Dec,2);
    Len_obj = size(Obj,2);
    Len_con = size(Con,2);
    for i = 1 : Len_obj
        [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
        [~,distinct2] = unique(round(Obj(:,i)*1e100)/1e100,'rows');
        distinct      = intersect(distinct1,distinct2);
        
        dmodel         = dacefit(Dec(distinct,:),Obj(distinct,i),'regpoly1','corrgauss',THETA_obj(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_obj{i}   = dmodel;
        THETA_obj(i,:) = dmodel.theta;
    end
    for i = 1 : Len_con
        [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
        [~,distinct2] = unique(round(Con(:,i)*1e100)/1e100,'rows');
        distinct      = intersect(distinct1,distinct2);
        
        dmodel         = dacefit(Dec(distinct,:),Con(distinct,i),'regpoly1','corrgauss',THETA_con(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_con{i}   = dmodel;
        THETA_con(i,:) = dmodel.theta;
    end
end

function [OffObj,Off_ObjMSE,OffCon,Off_ConMSE] = model_predict(Model_obj,Model_con,OffDec)
    N          = size(OffDec,1);
    Len_obj    = length(Model_obj);
    Len_con    = length(Model_con);
    OffObj     = zeros(N,Len_obj);
    OffCon     = zeros(N,Len_con);
    Off_ObjMSE = zeros(N,Len_obj);
    Off_ConMSE = zeros(N,Len_con);
    for i = 1 : N
        for j = 1 : Len_obj
            [OffObj(i,j),~,Off_ObjMSE(i,j)] = predictor(OffDec(i,:),Model_obj{j});
        end
        for j = 1 : Len_con
            [OffCon(i,j),~,Off_ConMSE(i,j)] = predictor(OffDec(i,:),Model_con{j});
        end
    end
    OffObj     = real(OffObj);
    OffCon     = real(OffCon);
    Off_ObjMSE = abs(real(Off_ObjMSE));
    Off_ConMSE = abs(real(Off_ConMSE));
end

function P = optimizaiton(Population,wmax,Model_obj,Model_con,Problem)
    global V 
    w      = 1;
    P.decs = Population.decs;
    while w <= wmax    
        OffDec   = OperatorGA(Problem,P.decs);
        P.decs   = [P.decs;OffDec];
        [P.objs,P.objmse,P.cons,P.conmse] = model_predict(Model_obj,Model_con,P.decs);
        index    = EnvironmentalSelection(P,V);
        P.decs   = P.decs(index,:);
        P.objs   = P.objs(index,:);
        P.cons   = P.cons(index,:);
        P.objmse = P.objmse(index,:);
        P.conmse = P.conmse(index,:);
        w        = w + 1;
    end
end