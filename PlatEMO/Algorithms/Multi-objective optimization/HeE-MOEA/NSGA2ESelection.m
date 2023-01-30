function Centers = NSGA2ESelection(tr_x, tr_y, models, str, Problem, Ke)
% NSGA-II based selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    nn=0;   
    while (nn==100||nn==0)
        population = rand(Problem.N,Problem.D); %N*D 
        population = population.*repmat(Problem.upper,Problem.N,1)+(1-population).*repmat(Problem.lower,Problem.N,1);%N*D   
        population(:,Problem.D+1:Problem.D+Problem.M)=Estimate(population(:,1:Problem.D), tr_x, tr_y, models, str, Problem.M);%The value of objective by model
              
        [~,FrontNo,CrowdDis] = EnvironmentalSelection(population,Problem.N,Problem.D,Problem.M);%N*(D+M)

        %% Optimization of NSGA-II
        for i=1:floor(Problem.maxFE/Problem.N)
            MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
            parent     = population(MatingPool,:);
            offspring  = OperatorGA(Problem,parent(:,1:Problem.D),{1,20,1,20});
            offspring(:,Problem.D+1:Problem.D+Problem.M)=Estimate(offspring(:,1:Problem.D), tr_x, tr_y, models, str, Problem.M);%The value of objective by model
            [population,FrontNo,CrowdDis] = EnvironmentalSelection([population;offspring],Problem.N,Problem.D,Problem.M);
        end
        [Centers, nn]=Kmean([population,FrontNo',CrowdDis'], tr_x, Problem.N, Problem.D, Ke);
    end
end

function objv_LCB=Estimate(x, tr_x, tr_y, models, str,m)
    alpha=2;
    TestSamNum=size(x,1);n=length(models);
    sum=zeros(TestSamNum, size(tr_y,2));
    for i=1:n
        clear y;
        stri=str(i,:);
        M=models(i).M;p=models(i).p;
        if strcmp(stri(1), 'FE')%PCA
            x_pca=x*M;
            if strcmp(stri(2), 'RBF1')%RBF
                y= sim(p,x_pca');y=y';
            elseif strcmp(stri(2), 'SVM')
                y=SVMtest(x_pca, p, TestSamNum,m);
            elseif strcmp(stri(2), 'RBF2')%RBF
                y=RBF2test(x_pca, p, TestSamNum);
            end

        elseif strcmp(stri(1), 'NONE')%NONE
            if strcmp(stri(2), 'RBF1')%RBF
                y= sim(p,x');y=y';
            elseif strcmp(stri(2), 'SVM')
                y=SVMtest(x, p, TestSamNum, m);
            elseif strcmp(stri(2), 'RBF2')%RBF
                y=RBF2test(x, p, TestSamNum);
            end

        elseif strcmp(stri(1), 'FS')%CSO
            x_cso=x(:,M);
            if strcmp(stri(2), 'RBF1')%RBF
                y= sim(p,x_cso');y=y';
            elseif strcmp(stri(2), 'SVM')
                y=SVMtest(x_cso, p, TestSamNum,m);
            elseif strcmp(stri(2), 'RBF2')%RBF
                y=RBF2test(x_cso, p, TestSamNum);
            end 
        end
        result(i).y=y;
        sum=sum+y;
    end
    me=sum/n;
    sum=zeros(TestSamNum, size(tr_y,2));
    for i=1:n
        y=result(i).y;
        sum=sum+(y-me).^2;
    end
    s2=sum/(n-1);
    %fmin=min(tr_y);
    %fmin_matrix=repmat(fmin,TestSamNum,1);
    objv_LCB=me-alpha*sqrt(s2);
end

function y=SVMtest(x, p, N, M)
    ps1=p{1};ps2=p{2};
    xtest0=mapminmax('apply',x',ps1);xtest0=xtest0';
%     py = zeros(N,M);
    for j=1:M
        for k=1:N
            py(k,j) = svmpredict(0,xtest0(k,:),p{j+2}, '-q');
        end
    end
    y=mapminmax('reverse',py',ps2);y=y';  
end

function y=RBF2test(x, p, N)
    Centers=p.Centers;Spreads=p.Spreads;
    W2=p.W2;B2=p.B2;
    TestDistance = dist(Centers,x');
    TestSpreadsMat = repmat(Spreads,1,N);
    TestHiddenUnitOut = radbas(TestDistance./TestSpreadsMat);
    y= (W2*TestHiddenUnitOut+repmat(B2,1,N))'; 
end

function [population,FrontNo,CrowdDis] = EnvironmentalSelection(population,N,D,M)
% The environmental selection of NSGA-II

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(population(:,D+1:D+M),[],N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(population(:,D+1:D+M),FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    population = population(Next,:);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

function [Centers, nn]=Kmean(chromosome, tr_x, N, D, Ke)
    % Delete the multiple solutions between chromosome and the train data
    Q=[];P=[];
    for i=1:N
        PandQ=[tr_x;Q];
        matrix=dist(chromosome(i,1:D),PandQ');%1*size(PandQ,1) array
        if isempty(find(matrix<=10^-06))
            Q=[Q;chromosome(i,1:D)];
            P=[P;chromosome(i,:)];
        end
    end
    if size(Q,1)<Ke
        Centers=[];nn=100;
    else
        index=randperm(size(Q,1));
        Centers=Q(index(1:Ke),1:D);n=1;si=size(Q,1);
        % Cluster by the distance of the decision space
        while n<100

            NumberInClusters = zeros(Ke,1); % Number of samples in each class ,default is 0
            IndexInClusters = zeros(Ke,N); % Index of samples in each class
            % Classify all samples by the least distance principle
            for i = 1:si
                AllDistance = dist(Centers,Q(i,1:D)');% Calculate the distance between the i-th solution and each clustering center
                [~,Pos] = min(AllDistance);   % Minimum distance,training input is the index of clustering center
                NumberInClusters(Pos) = NumberInClusters(Pos) + 1;
                IndexInClusters(Pos,NumberInClusters(Pos)) = i;% Stores the training indexes belonging to this class in turn
            end
            % Store the old clustering centers
            OldCenters = Centers;
            % Recalculate the clustering centers
            for i = 1:Ke
                Index = IndexInClusters(i,1:NumberInClusters(i));% Extract the training input index belonging to this class
                Centers(i,:) = mean(Q(Index,1:D),1);    %Take the average of each class as the new clustering center
            end
            % Judge whether the old and new clustering centers are consistent
            EqualNum = sum(sum(Centers==OldCenters));% Centers and OldCenters are subtracted from each other to sum up all corresponding bits
            if EqualNum ==D*Ke % The old and new clustering centers are consistent
                break,
            end
            n=n+1;
        end
        nn=n;fprintf('k-means clustering %d\n',n);
    end
end