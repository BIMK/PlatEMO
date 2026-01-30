classdef MiSACO < ALGORITHM
% <2022> <single> <real/integer> <expensive>
% Multi surrogate-assisted ant colony optimization
% Np ---  60 --- Size of parent population
% No --- 100 --- Size of offspring population

%------------------------------- Reference --------------------------------
% J. Liu, Y. Wang, G. Sun, and T. Pang. Multisurrogate-assisted ant colony
% optimization for expensive optimization problems with continuous and
% categorical variables. IEEE Transactions on Cybernetics, 2022, 52(11):
% 11348-11361.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiao Liu (email: jiao.liu@ntu.edu.sg)

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            [Np,No] = Algorithm.ParameterSet(60,100);
            PI      = UniformPoint(Np,Problem.D,'Latin');
            DB      = Problem.Evaluation(repmat(Problem.upper-Problem.lower,Np,1).*PI+repmat(Problem.lower,Np,1));
            P.decs  = DB.decs;
            P.decs(:,Problem.Dr+1:Problem.D) = ceil(P.decs(:,Problem.Dr+1:Problem.D));
            P.objs  = DB.objs;
            
            %% Optimization
            while Algorithm.NotTerminated(DB)
                %% Multisurrogate-Assisted Selection
                % Surrogate Construction
                Train_X = DB.decs;
                Train_X(:,Problem.Dr+1:Problem.D) = ceil(Train_X(:,Problem.Dr+1:Problem.D));
                Train_Y = DB.objs;
                upper_r = Problem.upper; upper_r = upper_r(:,1:Problem.Dr);
                lower_r = Problem.lower; lower_r = lower_r(:,1:Problem.Dr);
                RBF_glo = RBFCreate_Hamming(Train_X,Train_Y,Problem.Dc,upper_r,lower_r,'gaussian');
                RF_net  = fitrensemble(Train_X,Train_Y,'CategoricalPredictors',[Problem.Dr+1:Problem.D]);
                
                % Generate Offspring Population
                [~,rank_v] = sort(P.objs); [~,rank_v] = sort(rank_v);
                OffDec     = ACO_MV_generate(P.decs,Problem.Dr,Problem.Dc,Np,Problem.len_c,rank_v,No,upper_r,lower_r);
                
                % RBF-Based Selection,
                for i = 1 : No
                    RBF_evaluate(i) = RBFInterp_Hamming(OffDec(i,:),Problem.Dc, upper_r, lower_r, RBF_glo);
                end
                [~,rank] = sort(RBF_evaluate);
                DB       = [DB, Problem.Evaluation(OffDec(rank(1),:))];
                OffDec(rank(1),:) = [];
                
                % LSBT-Based Selection
                RF_evaluate = predict(RF_net,OffDec);
                [~,rank]    = sort(RF_evaluate);
                DB = [DB, Problem.Evaluation(OffDec(rank(1),:))];
                OffDec(rank(1),:) = [];
                
                % Random Selection
                DB = [DB, Problem.Evaluation(OffDec(randi(size(OffDec,1)),:))];
                
                %% Surrogate-Assisted Local Search
                [~,id_best] = min(DB.objs);
                if length(id_best) > 1
                    id_best = id_best(1,:);
                end
                best_x = DB(id_best).decs;
                best_r = best_x(:,1:Problem.Dr);
                best_c = ceil(best_x(:,Problem.Dr+1:Problem.D));
                Train_x    = DB.decs;
                idx_local  = dist_ham(best_c,ceil(Train_x(:,Problem.Dr+1:Problem.D)));
                best_Arc_r = Train_x(idx_local,1:Problem.Dr);
                best_Arc_f = DB(idx_local).objs;
                Local_flag = sum(idx_local)> 5*Problem.Dr;
                if Local_flag
                    idx_rank_best = dist_r(best_r,best_Arc_r);
                    local_up      = max(best_Arc_r(idx_rank_best,:));
                    local_dn      = min(best_Arc_r(idx_rank_best,:));
                    xr_Local      = Local_Search(best_Arc_r(idx_rank_best,:),best_Arc_f(idx_rank_best),local_up,local_dn,best_r);
                    DB            = [DB, Problem.Evaluation([xr_Local,best_c])];
                end
                
                % Population Update
                [~,rank] = sort(DB.objs);
                P.decs   = DB(rank(1:Np)).decs;
                P.decs(:,Problem.Dr+1:Problem.D) = ceil(P.decs(:,Problem.Dr+1:Problem.D));
                P.objs   = DB(rank(1:Np)).objs;
            end
        end
    end
end

function idx_sel = dist_ham(x,v)
    distt   = mean(x~=v,2);
    idx_sel = (distt == 0);
end

function out = dist_r(x0,x)
    dist = sqrt(sum((x0 - x).^2,2));
    [~,sort_rank] = sort(dist);
    [~,out_sort]  = sort(sort_rank);
    out = (out_sort<=30);
end