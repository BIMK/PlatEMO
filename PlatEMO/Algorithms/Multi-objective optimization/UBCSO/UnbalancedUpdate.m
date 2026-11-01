function NewCluster = UnbalancedUpdate(Problem,Cluster,numClusters,clusterSize,W,Z)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    Swarm            = [Cluster{:}];
    NewCluster       = cell(numClusters,1);
    LoserCluster     = cell(numClusters,1);
    Fitness          = CalFitness(Swarm.objs,Swarm.cons,0);
    N                = size(Swarm,2);
    [~,Association]  = min(acos(1-pdist2((Swarm.objs-repmat(Z,N,1)),W,'cosine')),[],2);
    SwarmNumberArray = zeros(1,numClusters);
    WinnerCluster    = cell(numClusters,1);
    for i = 1 : numClusters 
        ClusterSwarm        = Cluster{i};
        SwarmNumberArray(i) = length(ClusterSwarm);
        if i ~= numClusters
            ClusterFit = Fitness((i-1)*clusterSize+1:i*clusterSize);
        else
            ClusterFit = Fitness((i-1)*clusterSize+1:N);
        end
        if i ~= numClusters
            ClusterAssociate = Association((i-1)*clusterSize+1:i*clusterSize);
        else
            ClusterAssociate = Association((i-1)*clusterSize+1:N);
        end
        WorstGroup         = ClusterSwarm(ClusterAssociate~=i);
        c                  = length(WorstGroup);
        WLGroup            = ClusterSwarm(ClusterAssociate == i);
        WLGroupFitness     = ClusterFit(ClusterAssociate == i);
        WinnerGroup        = WLGroup(WLGroupFitness<1);
        WinnerGroupFitness = WLGroupFitness(WLGroupFitness<1);
        a                  = length(WinnerGroup);
        LoserGroup         = WLGroup(WLGroupFitness>=1);
        LoserGroupFitness  = WLGroupFitness(WLGroupFitness>=1);
        b                  = length(LoserGroup);

        if c==0 && b == SwarmNumberArray(i)
            [Loser,Winner] = LearnFromBest(Problem,LoserGroup,LoserGroupFitness);
        elseif c==0 && a == SwarmNumberArray(i)
            [Loser,Winner] = PairwiseCompetition(Problem,WinnerGroup,WinnerGroupFitness);
        elseif c==0 && (a+b) == SwarmNumberArray(i)
            [Loser,Winner] = LearnFromRandomWinners(Problem,WinnerGroup,LoserGroup);
        elseif c == SwarmNumberArray(i)
            Angle = acos(1-pdist2((WorstGroup.objs-repmat(Z,c,1)),W(i,:),'cosine'));
            [Loser,Winner] = LearnFromClosest(Problem,WorstGroup,Angle);
        else
            if (a+b) == 1
                [Loser,Winner] = LearnFromSingleWinner(Problem,WinnerGroup,LoserGroup,WorstGroup);
            else
                [Loser,Winner] = MixedCompetition(Problem,WinnerGroup,LoserGroup,WorstGroup,WinnerGroupFitness,LoserGroupFitness);
            end
        end
        WinnerCluster{i} = Winner;
        LoserCluster{i}  = Loser;
    end
    WinnerSwarm    = [WinnerCluster{:}];
    WinnerSwarmNum = length(WinnerSwarm);
    FitWS          = CalFitness(WinnerSwarm.objs,WinnerSwarm.cons,0);
    MatingPool     = TournamentSelection(2,WinnerSwarmNum,FitWS);
    OffspringWS    = OperatorGA(Problem,WinnerSwarm(MatingPool));
    if mod(WinnerSwarmNum,2) == 1
        OffspringWS = [OffspringWS OperatorGAhalf(Problem,[WinnerSwarm(randi(WinnerSwarmNum)) WinnerSwarm(randi(WinnerSwarmNum))])];
    end
    [~,AssociateWS] = min(acos(1-pdist2((OffspringWS.objs-repmat(Z,WinnerSwarmNum,1)),W,'cosine')),[],2);
    [~,Ind]         = sort(AssociateWS);
    OffspringWS     = OffspringWS(Ind);
    countWS         = 0;
    for i = 1 : numClusters
        No               = length(WinnerCluster{i});
        WinnerCluster{i} = OffspringWS(countWS+1:countWS+No);
        countWS          = countWS+No;
        NewCluster{i}    = [LoserCluster{i} WinnerCluster{i}];
    end
end

function [Loser,Winner] = LearnFromBest(Problem,Swarm,Fitness)
    MinFit       = min(Fitness);
    Index        = find(Fitness==MinFit);
    Index        = Index(randi(length(Index)));
    Winner       = Swarm(Index);
    a            = 1:length(Swarm);
    a(Index)     = [];
    Loser        = Swarm(a);
    NumWinner    = 1;
    LoserDec     = Loser.decs;
    WinnerDec    = Winner.decs;
    [NumLoser,D] = size(LoserDec);
    LoserVel     = Loser.adds(zeros(NumLoser,D));
    Winner.adds(zeros(NumWinner,D));
    r1           = repmat(rand(NumLoser,1),1,D);
    r2           = repmat(rand(NumLoser,1),1,D);
    OffVel       = r1.*LoserVel + r2.*(repmat(WinnerDec,NumLoser,1)-LoserDec);
    OffDec       = LoserDec + OffVel + r1.*(OffVel-LoserVel);
    OffDec       = PolynomialMutation(Problem,OffDec);
    Loser        = Problem.Evaluation(OffDec,OffVel);
end

function [Loser,Winner] = PairwiseCompetition(Problem,Swarm,Fitness)
    NumSwarm = size(Swarm,2);
    TWinner  = [];
    if mod(NumSwarm,2) == 1 
        MaxFit     = max(Fitness);
        Index      = find(Fitness==MaxFit);
        Index      = Index(randi(length(Index)));
        TWinner    = Swarm(Index);
        a          = 1:length(Swarm);
        a(Index)   = [];
        NewSwarm   = Swarm(a);
        NewFitness = Fitness(a);
    else
        NewSwarm   = Swarm;
        NewFitness = Fitness;
    end
    Rank           = randperm(length(NewSwarm),floor(length(NewSwarm)/2)*2);
    Loser          = Rank(1:end/2);
    Winner         = Rank(end/2+1:end);
    Change         = NewFitness(Loser) <= NewFitness(Winner);
    Temp           = Winner(Change);
    Winner(Change) = Loser(Change);
    Loser(Change)  = Temp;
    Winner         = NewSwarm(Winner);
    Loser          = NewSwarm(Loser);
    LoserDec       = Loser.decs;
    WinnerDec      = Winner.decs;
    [N,D]          = size(LoserDec);
    LoserVel       = Loser.adds(zeros(N,D));
    Winner.adds(zeros(N,D));
    r1             = repmat(rand(N,1),1,D);
    r2             = repmat(rand(N,1),1,D);
    OffVel         = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
    OffDec         = LoserDec + OffVel + r1.*(OffVel-LoserVel);

    OffDec = PolynomialMutation(Problem,OffDec);
    Loser  = Problem.Evaluation(OffDec,OffVel);
    Winner = [Winner TWinner];
end

function [Loser,Winner] = LearnFromRandomWinners(Problem,WinnerSwarm,LoserSwarm)
    TWinner        = WinnerSwarm;
    NumLoserSwarm  = length(LoserSwarm);
    Loser          = LoserSwarm;
    WinnerSwarmDec = WinnerSwarm.decs;
    LoserDec       = LoserSwarm.decs;
    WinnerDec      = zeros(NumLoserSwarm,size(LoserDec,2));
    for i = 1 : NumLoserSwarm
        No = randi(length(WinnerSwarm));
        WinnerDec(i,:) = WinnerSwarmDec(No,:);
    end
    [N,D]    = size(LoserDec);
    LoserVel = Loser.adds(zeros(N,D));
    r1       = repmat(rand(N,1),1,D);
    r2       = repmat(rand(N,1),1,D);
    OffVel   = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
    OffDec   = LoserDec + OffVel + r1.*(OffVel-LoserVel);

    %% Polynomial mutation
    OffDec = PolynomialMutation(Problem,OffDec);
    Loser  = Problem.Evaluation(OffDec,OffVel);
    Winner = TWinner;
end

function [Loser,Winner] = LearnFromClosest(Problem,Swarm,Fitness) 
    MinFit = min(Fitness);
    Index  = find(Fitness==MinFit);
    if ~isempty(Index)
        Index = Index(randi(length(Index)));
    else
        Index = 1;
    end
    Winner       = Swarm(Index);
    a            = 1 : length(Swarm);
    a(Index)     = [];
    Loser        = Swarm(a);
    NumWinner    = 1;
    LoserDec     = Loser.decs;
    WinnerDec    = Winner.decs;
    [NumLoser,D] = size(LoserDec);
    LoserVel     = Loser.adds(zeros(NumLoser,D));
    Winner.adds(zeros(NumWinner,D)); 
    r1           = repmat(rand(NumLoser,1),1,D);
    r2           = repmat(rand(NumLoser,1),1,D);
    OffVel       = r1.*LoserVel + r2.*(repmat(WinnerDec,NumLoser,1)-LoserDec);
    OffDec       = LoserDec + OffVel + r1.*(OffVel-LoserVel);
    OffDec       = PolynomialMutation(Problem,OffDec);
    Loser        = Problem.Evaluation(OffDec,OffVel);
end

function [Loser,Winner] = LearnFromSingleWinner(Problem,WinnerGroup,LoserGroup,WorstGroup) % a+b=1
    Winner       = [WinnerGroup LoserGroup];
    Loser        = WorstGroup;
    NumWinner    = 1;
    LoserDec     = Loser.decs;
    WinnerDec    = Winner.decs;
    [NumLoser,D] = size(LoserDec);
    LoserVel     = Loser.adds(zeros(NumLoser,D));
    Winner.adds(zeros(NumWinner,D)); 
    r1           = repmat(rand(NumLoser,1),1,D);
    r2           = repmat(rand(NumLoser,1),1,D);
    OffVel       = r1.*LoserVel + r2.*(repmat(WinnerDec,NumLoser,1)-LoserDec);
    OffDec       = LoserDec + OffVel + r1.*(OffVel-LoserVel);
    OffDec       = PolynomialMutation(Problem,OffDec);
    Loser        = Problem.Evaluation(OffDec,OffVel);
end

function [Loser,Winner] = MixedCompetition(Problem,WinnerGroup,LoserGroup,WorstGroup,WinnerGroupFitness,LoserGroupFitness) 
    [Loser1,~] = LearnFromRandomWinners(Problem,[WinnerGroup LoserGroup],WorstGroup);
    a          = length(WinnerGroup);
    b          = length(LoserGroup);
    c          = a + b;
    if b == c
        [Loser,Winner] = LearnFromBest(Problem,LoserGroup,LoserGroupFitness);
    elseif a == c
        [Loser,Winner] = PairwiseCompetition(Problem,WinnerGroup,WinnerGroupFitness);
    else
        [Loser,Winner] = LearnFromRandomWinners(Problem,WinnerGroup,LoserGroup);
    end
    Loser = [Loser Loser1];
end

function OffDec = PolynomialMutation(Problem,OffDec)
    [N,D]  = size(OffDec);
    Lower  = repmat(Problem.lower,N,1);
    Upper  = repmat(Problem.upper,N,1);
    disM   = 20;
    Site   = rand(N,D) < 1/D;
    mu     = rand(N,D);
    temp   = Site & mu<=0.5;
    OffDec = max(min(OffDec,Upper),Lower);
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end