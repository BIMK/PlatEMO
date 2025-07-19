classdef PEAplus < ALGORITHM
% <2024> <multi/many> <real/integer> <expensive> <constrained>
% Pareto-based Kriging-assisted constrained multiobjective evolutionary algorithm plus
% wmax --- 20 --- Generations of evolutionary search
% mu   ---  5 --- Number of selected candidates

%------------------------------- Reference --------------------------------
% Z. Zhang, Y. Wang, G. Sun, and K. Tang. Constrained probabilistic Pareto 
% dominance for expensive constrained multiobjective optimization problems. 
% IEEE Transactions on Evolutionary Computation, 2024.
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
            [wmax,mu] = Algorithm.ParameterSet(20,5);

            %% Initialization
            NI             = Problem.N;
            P              = UniformPoint(NI,Problem.D,'Latin');
            Database       = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            P              = Database; 
            THETA_obj      = 5.*ones(Problem.M,Problem.D);
            THETA_con      = 5.*ones(size(Database.cons,2),Problem.D);
            Model_obj      = cell(1,Problem.M);
            Model_con      = cell(1,size(Database.cons,2));
            sample_success = 1;

            %% Optimization
            while Algorithm.NotTerminated(Database)
                % Model Management
                if sample_success
                    [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(Database,THETA_obj,THETA_con);
                end
                % Evolutionary Search
                [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = Evo_Search(P,wmax,Model_obj,Model_con,Problem);
                % Candidate Seletion
                C = Candi_Select(PopDec,PopObj,PopCon,ObjMSE,ConMSE,Database,mu,Problem);
                % Population Update
                sample_success = 0;
                if isempty(C) == 0
                    C        = Problem.Evaluation(C);
                    Database = [Database,C];
                    index    = Pop_Update(Database.objs,Database.cons,NI);
                    P        = Database(index);
                    sample_success = 1;
                end
            end
        end
    end
end

function [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(Database,THETA_obj,THETA_con)
    Dec     = Database.decs;
    Obj     = Database.objs;
    Con     = Database.cons;
    Len_dec = size(Dec,2);
    Len_obj = size(Obj,2);
    Len_con = size(Con,2);
    for i = 1 : Len_obj
        [~,distinct1]  = unique(round(Dec*1e10)/1e10,'rows');
        [~,distinct2]  = unique(round(Obj(:,i)*1e10)/1e10,'rows');
        distinct       = intersect(distinct1,distinct2);
        dmodel         = dacefit(Dec(distinct,:),Obj(distinct,i),'regpoly1','corrgauss',THETA_obj(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_obj{i}   = dmodel;
        THETA_obj(i,:) = dmodel.theta;
    end
    for i = 1 : Len_con
        [~,distinct1]  = unique(round(Dec*1e10)/1e10,'rows');
        [~,distinct2]  = unique(round(Con(:,i)*1e10)/1e10,'rows');
        distinct       = intersect(distinct1,distinct2);
        dmodel         = dacefit(Dec(distinct,:),Con(distinct,i),'regpoly1','corrgauss',THETA_con(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_con{i}   = dmodel;
        THETA_con(i,:) = dmodel.theta;
    end
end