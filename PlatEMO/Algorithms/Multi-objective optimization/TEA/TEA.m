classdef TEA < ALGORITHM
% <2024> <multi/many> <real/integer> <expensive> <constrained/none>
% Two-phase evolutionary algorithm
% wmax --- 20 --- Generations of evolutionary search
% mu   ---  5 --- Number of selected candidates
    
%------------------------------- Reference --------------------------------
% Z. Zhang, Y. Wang, J. Liu, G. Sun, and K. Tang. A two-phase Kriging-
% assisted evolutionary algorithm for expensive constrained multiobjective
% optimization problems. IEEE Transactions on Systems, Man, and
% Cybernetics: Systems, 2024, 54(8): 4579-4591.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao_zhang0@163.com)
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            global phase NI Len_con Len_obj 
            [wmax,mu]      = Algorithm.ParameterSet(20,5,2);
            ct             = 0;
            ct_max         = 2;
            sample_success = 1;
            phase          = 1;

            %% Initialization
            NI        = Problem.N;
            P         = UniformPoint(NI,Problem.D,'Latin');
            Database  = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            P         = Database;
            Len_obj   = Problem.M;
            THETA_obj = 5.*ones(Len_obj,Problem.D);
            Model_obj = cell(1,Len_obj); 
            Len_con   = size(Database.cons,2);
            THETA_con = 5.*ones(Len_con,Problem.D);
            Model_con = cell(1,Len_con);

            %% Optimization
            while Algorithm.NotTerminated(Database)
                % Surrogate Construction
                if sample_success
                    [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(Database,THETA_obj,THETA_con);
                end
                % Evolutionary Search
                [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = Evo_Search(P,wmax,Model_obj,Model_con,Problem);
                % Candidate Seletion
                C = Candi_Select(PopDec,PopObj,PopCon,ObjMSE,ConMSE,Database,mu);
                % Phase Transition
                sample_success = 0;
                if isempty(C) == 0
                    C = Problem.Evaluation(C);
                    sample_success = 1;
                    [phase,ct]     = Phase_Trans(Database,C,ct,ct_max,phase);
                end
                % Population Reselection
                Database = [Database,C];
                index    = Pop_Reselect(Database.objs,Database.cons,NI);
                P        = Database(index);
            end
        end
    end
end

function [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(Database,THETA_obj,THETA_con)
    global Len_con Len_obj phase 
    Dec     = Database.decs;
    Obj     = Database.objs;
    Con     = Database.cons;
    Len_dec = size(Dec,2);
    Len_obj = size(Obj,2);
    Len_con = size(Con,2);
    for i = 1 : Len_obj
        [~,distinct1]  = unique(round(Dec*1e12)/1e12,'rows');
        [~,distinct2]  = unique(round(Obj(:,i)*1e12)/1e12,'rows');
        distinct       = intersect(distinct1,distinct2);
        dmodel         = dacefit(Dec(distinct,:),Obj(distinct,i),'regpoly1','corrgauss',THETA_obj(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_obj{i}   = dmodel;
        THETA_obj(i,:) = dmodel.theta;
    end
    if phase == 2
        for i = 1 : Len_con
            [~,distinct1]  = unique(round(Dec*1e12)/1e12,'rows');
            [~,distinct2]  = unique(round(Con(:,i)*1e12)/1e12,'rows');
            distinct       = intersect(distinct1,distinct2);
            dmodel         = dacefit(Dec(distinct,:),Con(distinct,i),'regpoly1','corrgauss',THETA_con(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
            Model_con{i}   = dmodel;
            THETA_con(i,:) = dmodel.theta;
        end
    else
        Model_con = [];
    end
end