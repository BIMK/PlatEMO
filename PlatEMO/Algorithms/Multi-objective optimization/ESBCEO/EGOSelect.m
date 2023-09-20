function PopDec = EGOSelect(Problem,Population,L1,L2,Ke,delta,nr)
% Selecting Points for Function Evaluation with the assistance of the
% Gaussian models

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    %% Fuzzy clustering the solutions
    [model,centers] = FCMmodel(Problem,Population,L1,L2);
    %% MOEA/D-DE optimization, where the popsize is set to N, the maximum evaluations is maxeva
	Z       = min(Population.objs,[],1);
    maxeva  = 20;
    eva=1;
    
	%% Generate the weight vectors
    [W, N]  = UniformPoint(Problem.N,Problem.M);
    T       = ceil(N/10);
    B       = pdist2(W,W);
    [~,B]   = sort(B,2);
    B       = B(:,1:T);
    duplicated = randi(length(Population),N,1);
    NewPop  = Population(duplicated);  % Sample N individuals from the initial population
    PopDec  = NewPop.decs; PopObj = NewPop.objs;
    gmin    = inf;

    [V0,N1] = UniformPoint(length(Population),Problem.M);
    V             = V0;
    alpha=2;
    fr=0.1;
    Population2=Population.decs;
    ObjValue2=Population.objs;
    while (eva<=maxeva)
    %%MOEAD/DE
        for i = 1 : N
            if rand < delta
                P = B(i,randperm(size(B,2)));
            else
                P = randperm(N);
            end
            OffDec = OperatorDE(Problem,PopDec(i,:),PopDec(P(1),:),PopDec(P(2),:));
            OffObj = Evaluate(Problem,OffDec,model,centers);
            Z = min(Z,OffObj);
            g_old = max(abs(PopObj(P,:) - repmat(Z,length(P),1)).*W(P,:),[],2);
            g_new = max(repmat(abs(OffObj-Z),length(P),1).*W(P,:),[],2);
            gmin = min([gmin,min(g_old),min(g_new)]);
            offindex = P(find(g_old>g_new,nr));
            if ~isempty(offindex)
                PopDec(offindex,:) = repmat(OffDec,length(offindex),1);
                PopObj(offindex,:) = repmat(OffObj,length(offindex),1);
            end
        end
        Pop1 = PopDec;
        Obj1 = PopObj;
        
        
       %% RVEA
        MatingPool = randi(length(Population2),1,N1);
        OffDec  = GABound(Population2(MatingPool',:),[Problem.upper;Problem.lower]);
        OffObj = Evaluate(Problem,OffDec,model,centers);
        PopDec=[Population2;OffDec];
        PopObj=[ObjValue2;OffObj];
        index = EnvironmentalSelection(ObjValue2,V,(Problem.FE/Problem.maxFE)^alpha);
        PopDec = PopDec(index,:);
        PopObj = PopObj(index,:);
        if ~mod(ceil(Problem.FE/N1),ceil(fr*Problem.maxFE/N1))
            V(1:N1,:) = ReferenceVectorAdaptation(PopObj,V0);
        end

        Pop2 = PopDec;
        Obj2 = PopObj;        
       
        if mod(eva,3)==0
           [FrontNo1,~]=NDSort(Obj1,inf);
           [FrontNo2,~]=NDSort(Obj2,inf);
            tag=1;
            PopDec=[Pop1(FrontNo1==1,:);Pop2(FrontNo2==1,:)];
            PopObj=[Obj1(FrontNo1==1,:);Obj2(FrontNo2==1,:)];
            if size(PopDec,1)>N
                d=zeros(size(PopDec,1),1);
                for i=1:size(PopDec,1)
                    for j=1:Problem.M
                        d(i)=d(i)+abs(PopObj(i,j))^(1/Problem.M);
                    end
                end
                dist=d.^(Problem.M);
                [dist,Index]=sort(dist);
                Dec=[];
                Obj=[];
                for j=1:N
                    Dec =[Dec;PopDec(Index==j,:)];
                    Obj =[Obj;PopObj(Index==j,:)];
                end
                PopDec=Dec;
                PopObj=Obj;

            else
                while size(PopDec,1)+size(Obj1(FrontNo1==tag+1,:),1)+size(Obj2(FrontNo2==tag+1,:),1)<N
                    tag=tag+1;
                    PopDec = [PopDec;Pop1(FrontNo1==tag,:);Pop2(FrontNo2==tag,:)];
                    PopObj = [PopObj;Obj1(FrontNo1==tag,:);Obj2(FrontNo2==tag,:)];
                end
                Dec=[Pop1(FrontNo1==tag+1,:);Pop2(FrontNo2==tag+1,:)];
                Obj=[Obj1(FrontNo1==tag+1,:);Obj2(FrontNo2==tag+1,:)];
                [Data,~]=ShanNonx([Dec Obj],[Dec Obj],size(Obj,2),10,10);
                [~,Index]=sort(Data,'descend');
                for j=1:N-size(PopDec,1)
                    PopDec =[PopDec;Dec(Index(j),:)];
                    PopObj =[PopObj;Obj(Index(j),:)];
                end
            end
        else
             PopDec=Pop1;
             PopObj=Obj1;  
        end
        Population2=PopDec;
        ObjValue2=PopObj;
        eva = eva +1;
    end
    
    %% Kmeans cluster the solutions into Ke clusters and select the solutions with the maximum EI in each cluster
    cindex  = kmeans(real(PopDec),Ke);
    Q = [];
    q=-0.5*cos(Problem.FE/Problem.maxFE*pi)+0.5;
    for i = 1 : Ke
        index = find(cindex == i); 
        temp = PopDec(index,:);
       
        [tempObj,~] = Evaluate(Problem,temp,model,centers);
        K = length(index);
        EI = zeros(K,1);
        [Data,~]=ShanNonx([temp tempObj],[temp tempObj],size(tempObj,2),10,10);
        for j = 1 : K
            EI(j) = EICal(tempObj(j,:),Data(j),1/size(tempObj,2),q);
        end
        [~,best] = max(EI);
        Q = [Q,index(best)];
    end
    PopDec = PopDec(Q,:);
end

function EI = EICal(Obj,s,p,q)
% Calculate the expected improvement

    M = size(Obj,2);
    d=0;
    for j=1:M
         d=d+abs(Obj(j))^p;
    end
    dist=d^(1/p);
    EI=(1-q)*dist-q*s;
end

function [y,x] = GPcal(lamda,mu,sig2)
% Calculate the mu (x) and sigma^2 (y) of the aggregation function

    tao = sqrt(lamda(1)^2*sig2(1) + lamda(2)^2*sig2(2));
    alpha = (mu(1)-mu(2))/tao;
    y = mu(1)*normcdf(alpha) + mu(2)*normcdf(-alpha) + tao*normpdf(alpha);
    x = (lamda(1)^2 + sig2(1))*normcdf(alpha) + ...
        (lamda(2)^2 + sig2(2))*normcdf(-alpha) + sum(lamda)*normpdf(alpha);
end

function [PopObj,MSE] = Evaluate(Problem,PopDec,model,centers)
% Predict the objective vector of the candidate solutions accodring to the
% Euclidean distance from each candidate solution to evaluated solutions
    
    D = pdist2(real(PopDec),real(centers));
    [~,index] = min(D,[],2);
    N = size(PopDec,1);
    PopObj = zeros(N,Problem.M);
    MSE = zeros(N,Problem.M);
    for i = 1 : N
        for j = 1 : Problem.M
            [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),model{index(i),j});
        end
    end
end

function [idx,dist] = nbselect(fitness,part,varargin)
    if varargin{1} == 'K'
        k = varargin{2};
        [idx,dist] = knnsearch(fitness(:,1:end),part(:,1:end),'Distance','euclidean','NSMethod','kdtree','K',k);  
    end
end