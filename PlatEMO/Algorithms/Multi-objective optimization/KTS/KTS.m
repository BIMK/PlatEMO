classdef KTS < ALGORITHM
% <multi/many> <real> <expensive> <constrained>
% Kriging-assisted evolutionary algorithm with two search modes
% tau --- 0.6 --- Threshold value
% phi --- 0.2 --- Threshold value
% mu  ---  20 --- Number of elite solution in A1

%------------------------------- Reference --------------------------------
% Z. Song, H. Wang, B. Xue, M. Zhang, and Y. Jin, Balancing objective
% optimization and constraint satisfaction in expensive constrained
% evolutionary multi-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2023.
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
            %% Parameter setting on KTA2
            phi1  = 0.1;
            wmax1 = 10;
            mu1   = 5;
            %% Parameter setting on KTS
            [tau,phi,mu] = Algorithm.ParameterSet(0.6,0.2,20);
            phi = -phi;
            
            %% Initialization
            p      = 1/Problem.M;
            CAsize = Problem.N;
            N      = Problem.N;
            P      = UniformPoint(N, Problem.D, 'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*P+repmat(Problem.lower,N,1));
            A1     = Population;
            
            % CA DA in KTA2
            CA = UpdateCA([],Population,CAsize);
            DA1 = UpdateDA(Population,[],Problem.N,p);
            
            % P1 P2 in KCCMO
            P1 = Population;
            P2 = Population;
            DA = DA1;
            
            % build surrogate models for objective functions
            Dec = [DA.decs;CA.decs;P2.decs];
            Obj = [DA.cons;CA.cons;P2.cons];
            [~,index] = unique( round(Dec*1e4)/1e4,'rows');
            Dec   =  Dec(index,:);
            Obj   = Obj(index,:);
            THETA = 5.*ones((Problem.M + size(DA.cons,2)),Problem.D);
            Model = cell(1,(Problem.M + size(DA.cons,2)));
            % build surrogate models for constraint functions
            for i = Problem.M+1: (Problem.M + size(DA.cons,2))
                j = i - Problem.M;
                dmodel     = dacefit(Dec,Obj(:,j),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                Model{i}   = dmodel;
                THETA(i,:) = dmodel.theta;
            end
            
            while Algorithm.NotTerminated(A1)
                Q          = Cal_Q(A1.objs); % calculate the convergence contribution of all evaluated solutions
                [Q, index] = sort(Q,'descend');
                CV         = sum(max(A1.cons,0),2); % calculate the CV of all evaluated solutions
                CV         = CV(index);
                coef = corrcoef(Q(end-mu:end), CV(end-mu:end)); % calculate the correlation coefficient
                r_coef = coef(2);
                % adaptive switching
                if r_coef < phi
                    search_mode = 1;   %  Constrained surrogate-assisted evolutionary search
                else
                    if r_coef <  tau
                        if rand < 0.5
                            search_mode = 0; %  Unconstrained surrogate-assisted evolutionary search
                        else
                            search_mode = 1;
                        end
                    else
                        search_mode = 0;
                    end
                end
                
                
                if search_mode == 0
                    DA = DA1;
                else
                    DA = P1;
                end
                
                % build surrogate models for objectives
                Dec = A1.decs;
                Obj = A1.objs;
                for i = 1:(Problem.M )%+ size(DA.cons,2))
                    dmodel     = dacefit(Dec,Obj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    THETA(i,:) = dmodel.theta;
                end
                
                if search_mode == 1
                    % build surrogate models for constraints
                    Dec = [DA.decs;CA.decs;P2.decs];
                    Obj = [DA.cons;CA.cons;P2.cons];
                    [~,index] = unique( roundn(Dec,-4) ,'rows');
                    Dec =  Dec(index,:);
                    Obj = Obj(index,:);
                    for i = Problem.M+1: (Problem.M + size(DA.cons,2))
                        j = i - Problem.M;
                        dmodel     = dacefit(Dec,Obj(:,j),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        Model{i}   = dmodel;
                        THETA(i,:) = dmodel.theta;
                    end
                end
                
                % Set the CCA and CDA as the current CA and DA
                CCA.obj = CA.objs; CCA.dec = CA.decs; CCA.con = CA.cons;  CCA.MSE = zeros(size(CCA.con,1),Problem.M+size(CCA.con,2));
                CP2.obj = P2.objs; CP2.dec = P2.decs; CP2.con = P2.cons; CP2.MSE = zeros(size(CP2.con,1),Problem.M+size(CP2.con,2));
                CDA.obj = DA.objs; CDA.dec = DA.decs; CDA.con = DA.cons; CDA.MSE = zeros(size(CDA.con,1),Problem.M+size(CDA.con,2));
                w = 1;
                while w <= wmax1
                    if search_mode == 0
                        % Solution generation in KTA2
                        [~,ParentCdec,~,ParentMdec] = MatingSelection_KTA2(CCA.obj,CCA.dec,CDA.obj,CDA.dec,Problem.N);
                        OffspringDec = [OperatorGA(Problem,ParentCdec,{1,20,0,0});OperatorGA(Problem,ParentMdec,{0,0,1,20})];
                    else
                        % Solution generation in KCCMO
                        Fitness1 = CalFitness(CP2.obj,CP2.con);
                        Fitness2 = CalFitness(CDA.obj);
                        %
                        MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                        MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                        %
                        Offspring1  = OperatorGA(Problem,CP2.dec(MatingPool1,:));
                        Offspring2  = OperatorGA(Problem,CDA.dec(MatingPool2,:));
                        OffspringDec = [Offspring1;Offspring2];
                    end
                    
                    Pop.dec = OffspringDec;
                    N       = size(Pop.dec,1);
                    Pop.obj = zeros(N,Problem.M);
                    Pop.con = zeros(N,size(DA.cons,2));
                    Pop.MSE = zeros(N,Problem.M+ size(DA.cons,2));
                    PopObj  = zeros(N,Problem.M+ size(DA.cons,2));
                    
                    for i = 1 : N
                        for j = 1 : (Problem.M+size(DA.cons,2))
                            [PopObj(i,j),~,Pop.MSE(i,j)] = predictor(Pop.dec(i,:),Model{j});
                        end
                    end
                    
                    Pop.obj = PopObj(:,1:Problem.M );
                    Pop.con = PopObj(:,Problem.M+1 :end );
                    
                    PopC = cat_struct(CCA,Pop);
                    CCA  = K_UpdateCA(PopC,CAsize);
                    PopD = cat_struct(CDA,Pop);
                    
                    if search_mode == 0
                        CDA = K_UpdateDA(PopD,Problem.N,p);
                    else
                        CDA = K_UpdateP(PopD,Problem.N,false);
                    end
                    
                    PopC1   = cat_struct(CP2,Pop);
                    [CP2,~] = K_UpdateP(PopC1,Problem.N,true);
                    w = w + 1;
                end
                
                % Adaptive sampling
                if search_mode == 0
                    % KTA2
                    % remove the same solution in all_population
                    [~,ia,~] = setxor(CCA.dec,A1.decs,'rows');
                    CCA = givevalue(CCA,ia);
                    [~,ia,~] = setxor(CDA.dec,A1.decs,'rows');
                    CDA = givevalue(CDA,ia);
                    
                    Offspring01 = Adaptive_sampling(CCA.obj,CDA.obj,CCA.dec,CDA.dec,CDA.MSE,DA,P2,mu1,p,phi1);
                else
                    % KCCMO
                    [CCA2,~]    = KCCMO_sampling(CP2,P2,mu1);
                    Offspring01 = CCA2.dec;
                end
                
                [~,index] = unique(Offspring01 ,'rows');
                PopNew = Offspring01(index,:);
                Offspring02 = PopNew;
                
                if ~isempty(Offspring02)
                    Offspring = Problem.Evaluation(Offspring02);
                    temp =  A1.decs;
                    for i = 1 : size(Offspring,2)
                        dist2 = pdist2(Offspring(i).decs,temp);
                        if min(dist2) > 1e-5
                            A1 = [A1,Offspring(i)];
                        end
                        temp = A1.decs;
                    end
                    
                    % data selection is the same as enviromental selection
                    CA     = UpdateCA(CA,Offspring,CAsize);
                    DA1    = UpdateDA(DA1,Offspring,Problem.N,p);
                    [P1,~] = Update_P([P1,Offspring],Problem.N,false);
                    [P2,~] = Update_P([P2,Offspring],Problem.N,true);
                end
            end
        end
    end
end

function value = Cal_Q(Obj)
    N = size(Obj,1);
    Obj = (Obj-repmat(min(Obj),N,1))./(repmat(max(Obj)-min(Obj),N,1));
    I = zeros(N);
    for i = 1 : N
        for j = 1 : N
            I(i,j) = max(Obj(i,:)-Obj(j,:));
        end
    end
    C = max(abs(I));
    F = sum(-exp(-I./repmat(C,N,1)/0.05)) + 1;
    value = 1./F;
end