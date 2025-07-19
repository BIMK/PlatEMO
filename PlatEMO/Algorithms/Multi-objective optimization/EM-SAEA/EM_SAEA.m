classdef EM_SAEA< ALGORITHM
% <multi/many> <real> <constrained> <expensive>
% Ensemble-based surrogate model-assisted evolutionary algorithm
% wmax   --- 20 --- Maximum number of generations before updating surrogate models
% lc_num ---  5 --- Number of local models

%------------------------------- Reference --------------------------------
% Y. Li, X. Feng, and H. Yu. Enhancing landscape approximation with
% ensemble-based surrogate model for expensive constrained multiobjective
% optimization. IEEE Transactions on Evolutionary Computation, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [wmax,lc_num] = Algorithm.ParameterSet(20,5); 

            %% Initialize weight vectors
            [V0,NI] = UniformPoint(Problem.N,Problem.M);
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            [ClW,~] = UniformPoint(lc_num,Problem.M);
            lc_num  = size(ClW,1);
            V       = V0;
            Vbb     = V0;
            T       = ceil(Problem.N/10);
            nr      = ceil(Problem.N/100);
            mu      = 5;
            mu1     = 1;
            kk      = 0.5;
            alpha   = 2;

            %% Detect the neighbours of each solution
            B     = pdist2(W,W);
            [~,B] = sort(B,2);
            B     = B(:,1:T);
            
            %% Initialize population by LHS
            PopDec     = UniformPoint(NI,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*PopDec+repmat(Problem.lower,NI,1));
            
            %% Initialize surrogate model settings
            numCon   = size(Population.cons,2);
            OTHETA   = 5.*ones(Problem.M,Problem.D);
            OModel   = cell(1,Problem.M);
            Model_lc = cell(lc_num,numCon);
            Model_gc = cell(1,numCon);
            THETA_lc = cell(lc_num);
            for i = 1 : lc_num
                THETA_lc{i} = 5.*ones(numCon,Problem.D);
            end
            THETA_gc = 5.*ones(numCon,Problem.D);  
            Z        = min(Population.objs,[],1);

            %% Calculate Angle Values
            angle      = acos(1-pdist2(W,W,'cosine'));
            temp_angle = angle;
            temp_angle(logical(eye(size(temp_angle)))) = inf;
            theta_min  = min(temp_angle');
            theta_min  = theta_min';
            theta      = theta_min.*0.5;
            TrainArc   = Population;
            A1         = TrainArc;
            stage      = 2;
            
            %% Optimization
            while Algorithm.NotTerminated(TrainArc)
                [~,index] = unique(TrainArc.decs,'rows');
                TrainArc  = TrainArc(index);
                TrainDec  = TrainArc.decs;
                TrainObj  = TrainArc.objs;
                TrainCon  = TrainArc.cons;
                [dsTrainDec,dsTrainObj] = dsmerge(TrainDec,TrainObj);
                for i = 1 : Problem.M
                    dmodel      = dacefit(dsTrainDec,dsTrainObj(:,i),'regpoly0','corrgauss',OTHETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    OModel{i}   = dmodel;
                    OTHETA(i,:) = dmodel.theta;
                end
                if stage == 1
                    %% Objective-oriented optimization stage(Code of RVMM from its authors)
                    % Refresh the model and generate promising solutions
                    A1Dec     = A1.decs;
                    A1Obj     = A1.objs;
                    [c,ia,ic] = unique(A1Obj,'rows','stable');
                    A1Obj     = A1Obj(ia,:);
                    A1Dec     = A1Dec(ia,:);
                    A1        = A1(ia);
	                
	                A1(sum(isnan(A1Obj),2)>0)      = [];
	                A1Dec(sum(isnan(A1Obj),2)>0,:) = [];
	                A1Obj(sum(isnan(A1Obj),2)>0,:) = [];
	                
	                A1(sum(isinf(A1Obj),2)>0)      = [];
	                A1Dec(sum(isinf(A1Obj),2)>0,:) = [];
	                A1Obj(sum(isinf(A1Obj),2)>0,:) = [];
                    
                    [c,ia,ic] = unique(A1Obj,'rows');
                    A1Obj     = A1Obj(ia,:);
                    A1Dec     = A1Dec(ia,:);
                    A1        = A1(ia);
                    [c,ib,ic] = unique(A1Dec,'rows');
                    A1Dec     = A1Dec(ib,:);
	                
                    V0         = Vbb;
                    PopDec1    = A1Dec;
                    w          = 1;
                    A2Obj      = TrainArc.objs;
                    [c,ia,ic]  = unique(A2Obj,'rows','stable');
                    A2Obj      = A2Obj(ia,:);
                    A2Obj      = A2Obj(NDSort(A2Obj,1)==1,:);
                    zmin       = min(A2Obj,[],1);
                    A2Obj_temp = A2Obj - repmat(min(zmin,[],1),size(A2Obj,1),1);
                    if size(A2Obj,1) >= 2
                        scale = (max(A2Obj,[],1)-min(A2Obj,[],1));
                        V_1   = V0.*(max(A2Obj,[],1)-min(A2Obj,[],1));
                    else
                        scale = ones(1,Problem.M);
                        V_1   = V0;
                    end
                    
                    Angle         = acos(1-pdist2(A2Obj_temp,V_1,'cosine'));
                    [~,associate] = min(Angle,[],2);
                    active        = unique(associate,'stable');
                    Va            = V_1(active,:);
                    NCluster      = min(5,size(Va,1));
                    [IDX,C]       = kmeans(V0(active,:),NCluster);
                    V1  = [];
                    ids = [];
                    for i = 1 : NCluster
                        EC      = find(IDX==i);
                        id      = EC(randperm(size(EC,1),1),1);
                        V1(i,:) = Va(id,:);
                        ids     = [ids;id];
                    end
                    if size(V1,1) < 5
                        notS = setdiff(1:size(V_1,1),active(ids));
                        V1   = [V1;V_1( notS(randperm(size(notS,2),5-size(V1,1))),:)];
                    end
                    WPopDec = [];
                    WPopObj = [];
                    WMSE    = [];
                    while w <= wmax
                        drawnow();
                        MatingPool = randi(size(PopDec1,1),1,Problem.N);
                        OffDec     = OperatorGA(Problem,PopDec1(MatingPool,:));
                        pop_candi  = [];
                        NP         = size(OffDec,1);
                        for ii = 1 : NP
                            if min(sqrt(sum((OffDec(ii,:) - [A1Dec;pop_candi]).^2,2)))>1E-6
                                pop_candi = [pop_candi;OffDec(ii,:)];
                            end
                        end
                        OffDec  = pop_candi;
                        PopDec1 = [PopDec1;OffDec];
                        [N,~]   = size(PopDec1);
                        PopObj1 = zeros(N,Problem.M);
                        MSE1    = zeros(N,Problem.M);
                        for i = 1 : N
                            for j = 1 : Problem.M
                                [PopObj1(i,j),~,MSE1(i,j)] = predictor(PopDec1(i,:),OModel{j});
                            end
                        end

                        MSE1 = max(MSE1,0);
                        S_   = sqrt(MSE1);
                        MSE1 = S_.*(MSE1<=1)+MSE1.*(MSE1>1);
		                
                        PopObj1_b = PopObj1;
                        MSE1_b    = MSE1;
                        PopObj1   = PopObj1+kk*MSE1;
                        
                        zmin = min([zmin;PopObj1],[],1);
                 
                        index = KEnvironmentalSelection(PopObj1,[V1;[]],(w/wmax)^alpha);
                        
                        PopDec1 = PopDec1(index,:);
                        PopObj1 = PopObj1(index,:);
                        MSE1    = MSE1(index,:);
                        
                        [~,ib]         = intersect(PopDec1,TrainArc.decs,'rows');
                        PopDec1_       = PopDec1;
                        PopDec1_(ib,:) = [];
                        if isempty(PopDec1_)
                            offobj1     = PopObj1_b(size(PopObj1_b,1)-size(OffDec,1)+1:end,:);
                            offmse1     = MSE1_b(size(PopObj1_b,1)-size(OffDec,1)+1:end,:);
                            [frontNo,~] = NDSort(offobj1,size(offobj1,1));
                            solId       = find(frontNo==1);
                            PopDec1     = [PopDec1; OffDec(solId,:)];
                            PopObj1     = [PopObj1; offobj1(solId,:)+kk*offmse1(solId,:)];
                            MSE1        = [MSE1; offmse1(solId,:)];
                        end
                        % Adapt referece vectors
                        if ~mod(w,ceil(wmax*0.1)) && size(unique(PopObj1,'rows'),1)>2
                            V1    = V1./scale;
                            scale = max(PopObj1,[],1)-min(PopObj1,[],1);
                            V1    = V1.*scale;
                        end
                        w = w + 1;
                    end
                    PopObj1 = PopObj1 - kk*MSE1;
                    
                    [c,ia,ic] = unique(PopObj1,'rows','stable');
                    if ~isempty(ia)
                        PopObj1 = PopObj1(ia,:);
                        PopDec1 = PopDec1(ia,:);
                        MSE1    = MSE1(ia,:);
                    end
                    
                    [~,ib]        = intersect(PopDec1,TrainArc.decs,'rows');
                    PopObj1(ib,:) = [];
                    PopDec1(ib,:) = [];
                    MSE1(ib,:)    = [];
                
                    w      = 1;
                    PopDec = A1Dec;
                    
                    while w <= wmax
                        drawnow();
                       
                        if size(PopDec,1) < Problem.N
                            MatingPool = randi(size(PopDec,1),1,Problem.N);
                            OffDec     = OperatorGA(Problem,PopDec(MatingPool,:));
                        else
                            OffDec = OperatorGA(Problem,PopDec);
                        end
                        
                        pop_candi = [];
                        NP        = size(OffDec,1);
                        for ii = 1 : NP
                            if min(sqrt(sum((OffDec(ii,:)-[A1Dec;pop_candi]).^2,2))) > 1E-6
                                pop_candi = [pop_candi;OffDec(ii,:)];
                            end
                        end
                        OffDec = pop_candi;
                        
                        PopDec = [PopDec;OffDec];
                        [N,~]  = size(PopDec);
                        PopObj = zeros(N,Problem.M);
                        MSE    = zeros(N,Problem.M);
                        for i = 1 : N
                            for j = 1 : Problem.M
                                [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),OModel{j});
                            end
                        end

                        MSE = max(MSE,0);
                        S_  = sqrt(MSE);
                        
                        MSE = S_.*(MSE<=1)+MSE.*(MSE>1);
		                
                        PopObj_b = PopObj;
                        PopDec_b = PopDec;
                        MSE_b    = MSE;
                        
                        PopObj = PopObj+kk*MSE;
                        
                        if w == 1
                            PopObj_ = PopObj(NDSort(PopObj,1)==1,:);
                            scale   = (max(PopObj_,[],1)-min(PopObj_,[],1));
                            scale(:,scale==0) = 10^(-6);
                            V = V0.*scale;
                            
                            Angle         = acos(1-pdist2(PopObj_, V,'cosine'));
                            [~,associate] = min(Angle,[],2);
                            active        = unique(associate,'stable');
                            Va            =  V(active,:);
                            if size(Va,1) < Problem.N
                                PopObj_temp = (PopObj_-min(PopObj_,[],1))./scale;
                                Vadd        = PopObj_temp(randperm(size(PopObj_temp,1),min(Problem.N-size(Va,1),size(PopObj_temp,1))),:);
                                V0          = [V0;Vadd];
                                V           =  V0.*scale;
                            end
                        end
                        
                        index   = KEnvironmentalSelection(PopObj,V,(w/wmax)^alpha);
                        PopDec  = PopDec(index,:);
                        PopObj  = PopObj(index,:);
                        MSE     = MSE(index,:);
                        [~,ib]  = intersect(PopDec,TrainArc.decs,'rows');
                        PopDec_ = PopDec;
                        PopDec_(ib,:) = [];
                        if isempty(PopDec_)
                            offobj      = PopObj_b(size(PopObj_b,1)-size(OffDec,1)+1:end,:);
                            offmse      = MSE_b(size(PopObj_b,1)-size(OffDec,1)+1:end,:);
                            [frontNo,~] = NDSort(offobj,size(offobj,1));
                            solId       = find(frontNo==1);
                            PopDec      = [PopDec; OffDec(solId,:)];
                            PopObj      = [PopObj;offobj(solId,:)+kk*offmse(solId,:)];
                            MSE         = [MSE; offmse(solId,:)];
                        end
                        
                        % Adapt referece vectors
                        if ~mod(w,ceil(wmax*0.1)) && size(unique(PopObj,'rows'),1)>2
                            V = V0.*(max(PopObj,[],1)-min(PopObj,[],1));
                        end
                        w       = w + 1;
                        WPopDec = [WPopDec;PopDec];
                        WPopObj = [WPopObj;PopObj];
                        WMSE    = [WMSE;MSE];            
                    end
	                 
                    PopObj = WPopObj;
                    PopDec = WPopDec;
                    MSE    = WMSE;
                    
                    PopObj = PopObj-kk*MSE;
                    
                    [c,ia,ic] = unique(PopObj,'rows','stable');
                    if ~isempty(ia)
                        PopObj = PopObj(ia,:);
                        PopDec = PopDec(ia,:);
                        MSE    = MSE(ia,:);
                    end
                    
                    [~,ib] = intersect(PopDec,TrainArc.decs,'rows');
                    
                    PopObj(ib,:) = [];
                    PopDec(ib,:) = [];
                    MSE(ib,:)    = [];
                    
                    NumVf  = [];
                    PopNew = KrigingSelect_RVMM(PopDec,PopObj,MSE,V,V0,NumVf,0.05*Problem.N,mu1,(w/wmax)^alpha,Problem.FE./Problem.maxFE,PopDec1,PopObj1,MSE1,V1,TrainArc,A2Obj);
                    
                    if ~isempty(PopNew)
                        [~,ib]       = intersect(PopNew,TrainArc.decs,'rows');
                        PopNew(ib,:) = [];
                        if ~isempty(PopNew)
                            New = Problem.Evaluation(PopNew);
                        else
                            New = [];
                        end
                    else
                        New = [];
                    end
                else
                    % Constraint-oriented optimization stage
                    ArcFZ      = min(TrainArc(all(TrainArc.cons<=0,2)).objs,[],1);
                    [PopDec,~] = CDPEnvironmentalSelection(TrainDec,TrainObj,TrainCon,Problem.N,W,ArcFZ);

                    % clustering
                    Zmin     = min(TrainObj,[],1);   
                    [N,~]    = size(TrainDec);
                    NormObj  = (TrainObj - Zmin);
                    Cluster1 = inf(N,lc_num);
                    Cluster2 = inf(N,lc_num);
                    temp     = N;
                    cluster_size = ceil(N/lc_num);
                    for i = 1 : lc_num
                        [~,ind] = sort(acos(1-pdist2(NormObj,ClW(i,:),'cosine')),"ascend");
                        Cluster1(ind(1:cluster_size),i) = i;
                    end
                    while temp > 0
                        for i = 1 : lc_num
                            [~,loc]         = min(acos(1-pdist2(NormObj,ClW(i,:),'cosine')));
                            Cluster2(loc,i) = i;
                            NormObj(loc,:)  = inf;
                            temp            = temp - 1;
                            if temp == 0
                                break;
                            end
                        end
                    end
                    Cluster = min(Cluster1,Cluster2);

                    % Train local models for constraints
                    for i = 1 : lc_num
                        for j = 1 : numCon
                            X_train_lc       = TrainDec(Cluster(:,i)==i,:); Y_train_lc = TrainCon(Cluster(:,i)==i,j);
                            [X_train_lc, Y_train_lc] = dsmerge(X_train_lc, Y_train_lc);
                            model_lc         = dacefit(X_train_lc,Y_train_lc,'regpoly0','corrgauss',THETA_lc{i}(j,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));                   
                            THETA_lc{i}(j,:) = model_lc.theta;
                            Model_lc{i,j}    = model_lc;
                        end
                    end

                    % Train global model for constraints
                    for i = 1 : numCon
                        [X_train_gc,Y_train_gc] = dsmerge(TrainDec,TrainCon(:,i));
                        model_gc      = dacefit(X_train_gc,Y_train_gc,'regpoly0','corrgauss',THETA_gc(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        THETA_gc(i,:) = model_gc.theta;
                        Model_gc{i}   = model_gc;
                    end

                    ArcZ   = min(TrainArc.objs,[],1);
                    ArcFZ  = min(TrainArc(all(TrainArc.cons<=0,2)).objs,[],1);
                    [N,~]  = size(PopDec);
                    PopObj = zeros(N,Problem.M);
                    MSE    = zeros(N,Problem.M);
                    for k = 1 : N
                        for j = 1 : Problem.M
                            [PopObj(k,j),~,MSE(k,j)] = predictor(PopDec(k,:),OModel{j});
                        end
                    end
                    PopCon = zeros(N,numCon);
                    CMSE   = zeros(N,numCon);
                    for k = 1 : N
                        for j = 1 : numCon
                            [PopCon(k,j),~,CMSE(k,j)] = predictor(PopDec(k,:),Model_gc{j});
                        end
                    end
                    PopCV = sum(max(0,PopCon),2);
                    w     = 0;
                    pf    = length(find(PopCV==0))/size(PopCV,1);
                    while w < wmax
                        for i = 1 : Problem.N
                            % Choose the parents
                            if rand < 0.9
                                P = B(i,randperm(size(B,2)));
                            else
                                P = randperm(Problem.N);
                            end
                            
                            % Generate an offspring
                            OffDec = OperatorGAhalf(Problem,PopDec(P(1:2),:));
                            for j = 1 : Problem.M
                                [OffObj(:,j),~,OffMSE(:,j)] = predictor(OffDec,OModel{j});
                            end
                            
                            % Update the ideal point
                            Zmin        = min(Zmin,OffObj);
                            Z           = min(Z,OffObj);
                            NormObj     = (OffObj - Zmin);
                            [~,loc]     = min(acos(1-pdist2(NormObj,ClW,'cosine')));
                            OffClus     = loc;
                            AnglePopObj = PopObj(P,:)-repmat(Z,length(P),1);
                            Angle0      = acos(1 - pdist2(real(AnglePopObj),W(P,:),'cosine'));
                            Angle       = diag(Angle0);
                            NewAngle    = acos(1-pdist2(real(OffObj-Z),W(P,:),'cosine'));
                            NewAngle    = NewAngle';
    
                            % Predict each constraint function value
                            for j = 1 : numCon
                                [OffCon_lc(:,j),~,OffMSE_lc(:,j)] = predictor(OffDec,Model_lc{OffClus,j});
                                [OffCon_gc(:,j),~,OffMSE_gc(:,j)] = predictor(OffDec,Model_gc{j});
                            end
                            if mean(OffMSE_gc) < mean(OffMSE_lc)
                                OffCon = OffCon_gc;
                            else
                                OffCon = OffCon_lc;
                            end
                        
                            % Calculate the constraint violation of offspring and P
                            CVO = sum(max(0,OffCon));
                            CVP = sum(max(0,PopCon(P,:)),2);
                            
                            % Calculate PBI value
                            normW   = sqrt(sum(W(P,:).^2,2));
                            normP   = sqrt(sum((PopObj(P,:)-repmat(Z,length(P),1)).^2,2));
                            normO   = sqrt(sum((OffObj-Z).^2,2));
                            CosineP = sum((PopObj(P,:)-repmat(Z,length(P),1)).*W(P,:),2)./normW./normP;
                            CosineO = sum(repmat(OffObj-Z,length(P),1).*W(P,:),2)./normW./normO;
                            g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                            g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
   
                            %% Vector-based Constrained Dominance Principle
                                if CVO+CVP == 0
                                    replace = find(g_old>=g_new,nr);
                                else
                                    if ~isempty(NewAngle<= theta(P) & Angle<= theta(P))
                                        replace = find((g_old>=g_new & CVP==CVO)| CVP>CVO,nr);
                                    else
                                        if rand<pf
                                            replace = find(g_old>=g_new,nr);
                                        else
                                            replace = find(CVP>CVO,nr);
                                        end
                                    end
                                end
                            if ~isempty(replace)
                                PopDec(P(replace),:) = OffDec;
                                PopObj(P(replace),:) = OffObj;
                                PopCon(P(replace),:) = OffCon;
                            end
                        end
                        w = w + 1;
                    end
                    % New solutions for re-evaluation
                    [NewDec,~] = ArchiveUpdate(PopDec,PopObj,PopCon,mu);
                    if ~isempty(NewDec)
                        New = Problem.Evaluation(NewDec);
                    else
                        NewDec = KrigingSelect(PopDec,PopObj,W,mu,(Problem.FE/Problem.maxFE)^2,PopCon);
                        New = Problem.Evaluation(NewDec);
                    end
                end
                TrainArc = [TrainArc,New];
                A1       = TrainArc;

                %% Stage switch
                if Problem.FE < ceil(NI+0.5*(Problem.maxFE-NI))
                    stage = 1;
                else
                    stage = 2;
                end
            end
        end
    end
end