function [Population,ObjBool] = WOF_Optimizer(Problem,Population,ConDim,DivDim,ConVars,DivVars,ObjBool)
% weighted optimization framework with NSGA-III

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Select diverse solution for WOF
    [XPrime,NumPrime] = SelectPrimeX(Population);
    PopSize = 10;

    % Set ideal point
    [Z,PopSize] = UniformPoint(PopSize,Problem.M);
    Zmin        = min(Population(all(Population.cons<=0,2)).objs,[],1);
    
    % Optimization
    for i = 1 : NumPrime
        % Create decision variable groups
        GroupIndex = VariableGrouping(Problem,XPrime(i),ConDim,DivDim,ConVars,DivVars);
        [WeightDecs,WOFPopulation] = InitWOFPopulation(Problem,GroupIndex,XPrime(i),PopSize,ConDim+DivDim);
        for j = 1 : 5
            FE = Problem.FE;
            MatingPool = TournamentSelection(2,2*PopSize,sum(max(0,WOFPopulation.cons),2));
            [OffspringWeight,Offspring] = WOFGAhalf(Problem,WeightDecs(MatingPool,:),GroupIndex,XPrime(i));
            Zmin = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
            [WOFPopulation,WeightDecs] = EnvironmentalSelection([WOFPopulation,Offspring],PopSize,Z,Zmin,[WeightDecs;OffspringWeight]);
            ObjBool(FE+1:Problem.FE,:) = Offspring.objs;
        end
        Population = WOF_EnvironmentalSelection([Population,WOFPopulation],Problem.N);
    end
end

function GroupIndex = VariableGrouping(Problem,XPrime,ConDim,DivDim,ConVars,DivVars)
    
    % Convergence variable: Order By Value Grouping
    NumCon     = numel(ConVars);
    NumPerCon  = floor(NumCon/ConDim);
    [~,Index]  = sort(XPrime.dec(ConVars),'ascend');
    GroupIndex = ones(1,Problem.D);
    for i = 1 : ConDim-1
        GroupIndex(ConVars(Index(((i-1)*NumPerCon)+1:i*NumPerCon))) = i;
    end
    GroupIndex(ConVars(Index(((ConDim-1)*NumPerCon)+1:end))) = ConDim;

    % Diversity variable: Order By Value Grouping
    NumDiv    = numel(DivVars);
    NumPerCon = floor(NumDiv/DivDim);
    [~,Index] = sort(XPrime.dec(DivVars),'ascend');
    for i = 1 : DivDim-1
        GroupIndex(DivVars(Index(((i-1)*NumPerCon)+1:i*NumPerCon))) = i+ConDim;
    end
    GroupIndex(DivVars(Index(((DivDim-1)*NumPerCon)+1:end))) = ConDim+DivDim;
end

function [XPrime,NumPrime] = SelectPrimeX(Population)
    % m+1 by reference lines
    M        = size(Population.objs,2);
    NumPrime = M + 1;
    Vectors  = [eye(M);ones(1,M)];
    XPrime   = [];
    for i = 1 : NumPrime
        Angles    = pdist2(Vectors,Population.objs,'cosine');
        [~,Index] = min(Angles);
        XPrime    = [XPrime,Population(Index)];
    end
end

function [WeightDecs,WOFPopulation] = InitWOFPopulation(Problem,GroupIndex,XPrime,PopSize,DimReduce)
    WeightDecs    = rand(PopSize,DimReduce).*2.0;
    PopDecs       = ReverseDecs(Problem,XPrime,WeightDecs(:,GroupIndex));
    WOFPopulation = Problem.Evaluation(PopDecs);
end

function PopDecs = ReverseDecs(Problem,XPrime,WeightDecs)
% The interval-intersection method

    [PopSize,~] = size(WeightDecs);
    XPrime      = repmat(XPrime.decs,PopSize,1);
    LowerBound  = repmat(Problem.lower,PopSize,1);
    UpperBound  = repmat(Problem.upper,PopSize,1);

    Interval = XPrime - LowerBound;
    PopDecs  = LowerBound + WeightDecs.*Interval;
    Interval = UpperBound - XPrime;
    PopDecs(WeightDecs>1) = XPrime(WeightDecs>1) + (WeightDecs(WeightDecs>1)-1).*Interval(WeightDecs>1);
    PopDecs  = min(max(PopDecs,LowerBound),UpperBound);
end

function [OffspringWeight,Offspring] = WOFGAhalf(Problem,Parent,GroupIndex,XPrime)
% GA for WOF

    %% Parameter setting
    [proC,disC,proM,disM] = deal(1,20,1,20);
    Parent1 = Parent(1:floor(end/2),:);    
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);

    %% Genetic operators for real encoding
    % Simulated binary crossover
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;

    % Polynomial mutation
    Lower = repmat(zeros(1,D),N,1);
    Upper = repmat(2*ones(1,D),N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    
    OffspringWeight = Offspring;
    OffDecs         = ReverseDecs(Problem,XPrime,OffspringWeight(:,GroupIndex));
    Offspring       = Problem.Evaluation(OffDecs);
end

function [Population,WeightDecs] = EnvironmentalSelection(Population,N,Z,Zmin,WeightDecs)
% The environmental selection of NSGA-III

    if nargin > 4
        tag = true;
    else
        tag = false;
    end

    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    Choose = LastSelection(Population(Next).objs,Population(Last).objs,N-sum(Next),Z,Zmin);
    Next(Last(Choose)) = true;
    % Population for next generation
    Population = Population(Next);

    if tag
        WeightDecs = WeightDecs(Next,:);
    else
        WeightDecs = [];
    end
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end

function [Population,FrontNo,CrowdDis] = WOF_EnvironmentalSelection(Population,N)
% The environmental selection of NSGA-II that is used inside WOF.
    
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = false(1,length(FrontNo));
    Next(FrontNo<MaxFNo) = true;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Popsize = min(N,size(Population,2));
    Next(Last(Rank(1:Popsize-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end