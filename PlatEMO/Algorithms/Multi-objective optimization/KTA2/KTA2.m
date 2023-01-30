classdef KTA2 < ALGORITHM
% <multi/many> <real/integer> <expensive>
% Kriging-assisted Two_Arch2
% tau  --- 0.75 --- Proportion of one type noninfluential points in training data
% phi  ---  0.1 --- Number of randomly selected individuals
% wmax ---   10 --- Number of generations before updating CA and DA 
% mu   ---    5 --- Number of re-evaluated solutions at each generation

%------------------------------- Reference --------------------------------
% Z. Song, H. Wang, C. He and Y. Jin, A Kriging-assisted two-archive
% evolutionary algorithm for expensive many-objective optimization, IEEE
% Transactions on Evolutionary Computation, 2021, 25(6): 1013-1027.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhenshou Song
% Email:zssong@stu.xidian.edu.cn

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [tau,phi,wmax,mu] = Algorithm.ParameterSet(0.75,0.1,10,5);
            
            %% Initialization 
            p = 1/Problem.M;
            CAsize = Problem.N;
            N = Problem.N;
            P = UniformPoint(N, Problem.D, 'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*P+repmat(Problem.lower,N,1));
            All_Population = Population;
            Ho_Population = All_Population;
            CA = UpdateCA([],Population,CAsize);
            DA = Population;
            THETA_S = 5.*ones(Problem.M,Problem.D);
            THETA_IS =  5.*ones(2,Problem.M,Problem.D);
            Model_sensitive = cell(1,Problem.M);
            Model_insensitive = cell(2,Problem.M);
            %% Optimization
            while Algorithm.NotTerminated(All_Population)
                %***** Building influential point-insensitive model********
                % build sensitive model
                Dec = All_Population.decs;
                Obj = All_Population.objs;
                for i = 1:Problem.M
                    dmodel     = dacefit(Dec,Obj(:,i),'regpoly0','corrgauss',THETA_S(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model_sensitive{i}   = dmodel;
                    THETA_S(i,:) = dmodel.theta;
                end
                % build insensitive models 
                Centers = zeros(Problem.M,2);
                for i = 1 : Problem.M
                    [~,N1] = sort(Obj(:,i));
                    num = ceil(length(All_Population).*tau);
                    mean_index{1} = N1(1:num);
                    mean_index{2} = N1(end-num:end);
                    for j = 1:2
                        Centers(i,j) = mean(Obj(mean_index{j},i));  % lambda and miu
                    end
                    for j = 1 : 2
                        train_X = Dec(mean_index{j},:);
                        train_Y = Obj(mean_index{j},i);
                        dmodel  = dacefit(train_X,train_Y,'regpoly0','corrgauss',THETA_IS(j,i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        Model_insensitive{j,i} = dmodel;
                        THETA_IS(j,i,:)        = dmodel.theta;
                    end
                end
                % Set the CCA and CDA as the current CA and DA
                CAobj = CA.objs; CAdec = CA.decs;
                DAobj = DA.objs; DAdec = DA.decs;
                w = 1;
                while w <= wmax   % this part is same as Two_Arch2 
                    [~,ParentCdec,~,ParentMdec] = MatingSelection_KTA2(CAobj,CAdec,DAobj,DAdec,Problem.N);
                    OffspringDec = [OperatorGA(Problem,ParentCdec,{1,20,0,0});OperatorGA(Problem,ParentMdec,{0,0,1,20})];
                    PopDec = [DAdec;CAdec;OffspringDec];
                    N      = size(PopDec,1);
                    PopObj = zeros(N,Problem.M);
                    MSE    = zeros(N,Problem.M);
                    %****** Using influential point-insensitive model *****
                    for i = 1:N
                        for j = 1:Problem.M
                            [PopObj(i,j),~,~] = predictor(PopDec(i,:),Model_sensitive{j});
                            if abs(PopObj(i,j)- Centers(j,1)) <= abs(PopObj(i,j)- Centers(j,2))
                                model = Model_insensitive{1,j};
                            else
                                model = Model_insensitive{2,j};
                            end
                            [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),model);
                        end
                    end
                    [CAobj,CAdec,~] = K_UpdateCA(PopObj,PopDec,MSE,CAsize);
                    [DAobj,DAdec,DAvar] = K_UpdateDA(PopObj,PopDec,MSE,Problem.N,p);
                    w = w + 1;
                end
                
                % Adaptive sampling 
                Offspring01 = Adaptive_sampling(CAobj,DAobj,CAdec,DAdec,DAvar,DA,mu,p,phi);
                
                [~,index] = unique(Offspring01 ,'rows');
                PopNew = Offspring01(index,:);
                Offspring02 = [];
                for i = 1:size(PopNew,1)
                    dist2 = pdist2(real( PopNew(i,:)),real(All_Population.decs));
                    if min(dist2) > 1e-5
                        Offspring02 = [Offspring02;PopNew(i,:)];
                    end
                end
                if ~isempty(Offspring02)
                    Offspring = Problem.Evaluation(Offspring02);

                    temp =  All_Population.decs;
                    for i = 1:size(Offspring,2)
                        dist2 = pdist2(real(Offspring(i).decs),real(temp));
                        if min(dist2) > 1e-5
                            All_Population = [All_Population,Offspring(i)];
                        end
                        temp = All_Population.decs;
                    end
                    CA = UpdateCA(CA,Offspring,CAsize);
                    DA = UpdateDA(DA,Offspring,Problem.N,p);
                    Ho_Population = [Ho_Population,Offspring];
                end
            end
        end
    end
end