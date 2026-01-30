classdef SAMOEATL2M < ALGORITHM
% <2025> <multi/many> <real/integer> <expensive>
% Surrogate-assisted multiobjective evolutionary algorithm based on two-level model management
% G     ---  20 --- Number of generations before updating suoorgate models
% KE    ---   5 --- The maximum number of infill solutions at each generation
% alpha --- 0.4 --- The threshold for the accuracy rate

%------------------------------- Reference --------------------------------
% Y. Liu, J. Ding, Q. Li, F. Li, and J. Liu. A two-level model
% management-based surrogate-assisted evolutionary algorithm for
% medium-scale expensive multiobjective optimization. IEEE Transactions on
% Systems, Man, and Cybernetics: Systems, 2025, 55(11): 8166-8180.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yuanchao Liu (email:liuyuanchao@ise.neu.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [G,KE,alpha] = Algorithm.ParameterSet(20,5,0.4);
            stage = 1;

            %% Generate the reference points and population
            [~,Problem.N] = UniformPoint(Problem.N,Problem.M);
            NI  = 100;
            NP  = NI;
            NI  = 11*Problem.D-1;
            P   = UniformPoint(NI, Problem.D, 'Latin');
            A2  = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            Num = randi([1, NI], 1, NP);
            A1  = A2(Num);

            %% Optimization
            while Algorithm.NotTerminated(A2)
                nobj   = Problem.M;
                yyy    = A2.objs;
                x      = A2.decs;
                yy     = (yyy-min(yyy))./(max(yyy)-min(yyy));
                PopObj = A1.objs;
                PopObj = (PopObj-min(yyy))./(max(yyy)-min(yyy));
                for i = 1 : nobj
                    dmodel   = TrainModel2(x,yy(:,i),Problem.M,Problem.D);
                    Model{i} = dmodel;
                end
                Combine = A2;
                fmin    = repmat(min(A2.objs,[],1),length(Combine),1);
                fmax    = repmat(max(A2.objs,[],1),length(Combine),1);
                PopObjA = (Combine.objs-fmin)./(fmax-fmin);
                [N,~]   = size(PopObjA);
                % Calculate the shifted distance between each two solutions
                sde = inf(N);
                for k = 1 : N
                    SPopObj = max(PopObjA,repmat(PopObjA(k,:),N,1));
                    for j = [1:k-1,k+1:N]
                        sde(k,j) = norm(PopObjA(k,:)-SPopObj(j,:));
                    end
                end
                SDE    = min(sde,[],2);
                Cd     = (SDE-min(SDE))./(max(SDE)-min(SDE));
                x      = A2.decs;
                dmodel = TrainModel2(x,Cd,Problem.M,Problem.D);
                ModelC = dmodel;
                % RBF-Assisted Evolutionary Multi-Objective Search
                PopDec = A1.decs;
                w      = 1;
                while w <= G
                    drawnow();
                    OffDec = OperatorGA(Problem,PopDec);
                    PopDec = [PopDec;OffDec];
                    OffObj = [];
                    for i = 1 : nobj
                        TestSamNum  = size(OffDec,1);
                        OffObj(:,i) = RBF(OffDec, Model{i}, TestSamNum);
                    end
                    PopObj = [PopObj;OffObj];
                    [FrontNo,MaxFNo] = NDSort(PopObj,size(OffDec,1));
                    Next       = FrontNo < MaxFNo;
                    TestSamNum = size(PopDec,1);
                    PopD       = RBF(PopDec, ModelC, TestSamNum);
                    Last       = find(FrontNo==MaxFNo);
                    [~,Rank]   = sort(PopD(Last),'descend');
                    Next(Last(Rank(1:size(OffDec,1)-sum(Next)))) = true;
                    PopObj = PopObj(Next,:);
                    PopDec = PopDec(Next,:);
                    w      = w + 1;
                end
                PopDec1 = PopDec;
                % Two-Level Model Management Strategy
                Infill_dataC  = [];
                dist2         = pdist2(real(PopDec1),real(A2.decs));
                DF            = find(min(dist2,[],2)<1e-5);
                PopDec1(DF,:) = [];
                PopObj(DF,:)  = [];
                if ~isempty(PopDec1)
                    if stage == 2
                        NumC = round(KE*(1-alpha));
                        if NumC < 1
                            NumC = 1;
                        end
                        [FrontNo,MaxFNo] = NDSort(PopObj,1);
                        Last = find(FrontNo==MaxFNo);
                        s    = [];
                        Pdec = A2.decs;
                        Pobj = A2.objs;
                        for h = 1 : Problem.M
                            [~, s(:,h)] = idw_prediction_and_uncertainty(PopDec1, Pdec, Pobj(:,h));
                        end
                        PopD         = sum(s,2);
                        [~,Rank]     = sort(PopD(Last),'descend');
                        Infill_dataC = [];
                        if size(PopDec1,1) > NumC
                            h = 1;
                            while size(Infill_dataC,1) < NumC
                                if h == 1
                                    Infill_dataC = PopDec1(Last(Rank(h)),:);
                                else
                                    Infill_TdataC = PopDec1(Last(Rank(h)),:);
                                    dist2 = pdist2(real(Infill_TdataC),real(Infill_dataC));
                                    if min(dist2) > 1e-5
                                        Infill_dataC = [Infill_dataC;Infill_TdataC];
                                    end
                                end
                                h = h +1;
                            end
                        else
                            Infill_dataC = PopDec1;
                        end
                    end
                    if stage ==  1
                        Num = KE;
                    else
                        Num = KE - NumC;
                    end
                    [FrontNo,MaxFNo] = NDSort(PopObj,1);
                    Last = find(FrontNo==MaxFNo);
                    if length(Last) <= Num
                        Infill_data = PopDec1(Last,:);
                    else
                        Next    = FrontNo < MaxFNo;
                        Combine = A2.objs;
                        Combine = [Combine;PopObj];
                        fmin    = repmat(min(Combine,[],1),size(Combine,1),1);
                        fmax    = repmat(max(Combine,[],1),size(Combine,1),1);
                        PopObjA = (Combine-fmin)./(fmax-fmin);
                        [N,~]   = size(PopObjA);
                        % Calculate the shifted distance between each two solutions
                        sde     = inf(N);
                        SPopObj = [];
                        for k = 1 : N
                            SPopObj = max(PopObjA,repmat(PopObjA(k,:),N,1));
                            for j = [1:k-1,k+1:N]
                                sde(k,j) = norm(PopObjA(k,:)-SPopObj(j,:));
                            end
                        end
                        SDE         = min(sde,[],2);
                        Cd          = (SDE-min(SDE))./(max(SDE)-min(SDE));
                        PopD        = Cd(length(A2)+1:end);
                        [~,Rank]    = sort(PopD(Last),'descend');
                        Infill_data = [];
                        if size(PopDec1,1) > Num
                            h = 1;
                            while size(Infill_data,1) < Num
                                if h == 1
                                    Infill_data = PopDec1(Last(Rank(h)),:);
                                else
                                    Infill_Tdata = PopDec1(Last(Rank(h)),:);
                                    dist2 = pdist2(real(Infill_Tdata),real(Infill_data));
                                    if min(dist2) > 1e-5
                                        Infill_data = [Infill_data;Infill_Tdata];
                                    end
                                end
                                h = h +1;
                            end
                        else
                            Infill_data = PopDec1;
                        end
                    end
                    if stage ==  2
                        Infill_data = [Infill_data;Infill_dataC];
                    end
                    % Update population
                    if ~isempty(Infill_data)
                        New     = Problem.Evaluation(Infill_data);
                        A2      = [A2,New];
                        NextA1  = [A1,New];
                        [FrontNoNA1,MaxFNo] = NDSort(NextA1.objs,NextA1.cons,NP);
                        Next    = FrontNoNA1 < MaxFNo;
                        Combine = NextA1.objs;
                        fmin    = repmat(min(Combine,[],1),size(Combine,1),1);
                        fmax    = repmat(max(Combine,[],1),size(Combine,1),1);
                        PopObjA = (Combine-fmin)./(fmax-fmin);
                        [N,~]   = size(PopObjA);
                        % Calculate the shifted distance between each two solutions
                        sde = inf(N);
                        SPopObj = [];
                        for k = 1 : N
                            SPopObj = max(PopObjA,repmat(PopObjA(k,:),N,1));
                            for j = [1:k-1,k+1:N]
                                sde(k,j) = norm(PopObjA(k,:)-SPopObj(j,:));
                            end
                        end
                        SDE      = min(sde,[],2);
                        Cd       = (SDE-min(SDE))./(max(SDE)-min(SDE));
                        Last     = find(FrontNoNA1==MaxFNo);
                        [~,Rank] = sort(Cd(Last),'descend');
                        Next(Last(Rank(1:NP-sum(Next)))) = true; %%Use the shifted distance to select the offspring
                        A1       = NextA1(Next);
                    end
                    % ARI
                    if ~isempty(Infill_data)
                        rp        = sum(Next(NP+1:end))/size(Infill_data,1);
                        FrontNoA1 = NDSort(A1.objs,A1.cons,NP);
                        rn        = length(find(FrontNoA1(NP - sum(Next(NP+1:end))+1:end)==1))/sum(Next(NP+1:end));
                        ARI       = min(rp, rn);
                        if ARI < alpha
                            stage = 2;
                        else
                            stage = 1;
                        end
                    else
                        stage = 2;
                    end
                else
                    NewPop = UniformPoint(1, Problem.D, 'Latin');
                    New    = Problem.Evaluation(NewPop);
                    A2     = [A2,New];
                end
            end
        end
    end
end

function y = RBF(x, p, N)
    Centers = p.Centers;
    Spreads = p.Spreads;
    W2      = p.W2;
    B2      = p.B2;
    TestDistance      = dist(Centers,x');
    TestSpreadsMat    = repmat(Spreads,1,N);
    TestHiddenUnitOut = radbas(TestDistance./TestSpreadsMat);
    y = (W2*TestHiddenUnitOut+repmat(B2,1,N))';
end