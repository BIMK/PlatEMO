classdef UBCSO < ALGORITHM
% <2026> <multi> <real> <constrained>
% Unbalanced and balanced competition strategy-assisted dual-swarm optimizer

%------------------------------- Reference --------------------------------
% Y. Wang, C. Hu, X. Wu, X. Yan, W. Gong, and L. Gao. Unbalanced and
% balanced competition strategy-assisted dual-swarm optimizer for
% constrained multi-objective optimization. IEEE Transactions on Automation
% Science and Engineering, 2026, 23: 7572-7586.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Generate the sampling points and random population
            if Problem.M == 3 && Problem.N ~= 90
                Problem.N = 90;
            end
            clusterSize     = 2;
            [W,~]           = UniformPoint(floor(Problem.N/clusterSize),Problem.M);
            numClusters     = size(W,1);
            Swarm2          = Problem.Initialization();
            Swarm1          = Problem.Initialization();
            Z               = min(Swarm2.objs,[],1);
            changeThreshold = 1e-1;
            stage           = 1; 
            maxChange       = 1;
            IdealHistory    = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            NadirHistory    = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            Archive         = ArchiveUpdate(Swarm2,Problem.N);
            generation      = ceil(Problem.FE/(Problem.N*2));
            G               = Problem.maxFE/(Problem.N*2);
            historyLength   = 50;
            threshold       = 0.5*G;
            transitionGen   = 0;
            Fitness1        = CalFitness(Swarm1.objs,Swarm1.cons,inf);
            mergeCount      = 1;
            mergeInterval   = floor(G/4);
            Swarm2          = AssociateSwarm(Swarm2,W,Z);
            Cluster         = ClusterDivision(Swarm2,numClusters,clusterSize);

            %% Optimization
            while Algorithm.NotTerminated(Swarm2)
                IdealHistory(generation,:) = Z;
                swarm                      = [Swarm1.decs,Swarm1.objs];
                NadirHistory(generation,:) = max(swarm(:,Problem.D + 1 : Problem.D + Problem.M),[],1);
                if generation >= threshold
                    maxChange = MaximumChange(IdealHistory,NadirHistory,generation,historyLength);
                end
                if maxChange <= changeThreshold && stage == 1
                    stage         = -1;
                    Fitness1      = CalFitness(Swarm1.objs,Swarm1.cons,0);
                    remainingGen  = G-generation;
                    mergeInterval = floor(remainingGen/4);
                    transitionGen = generation;
                    Swarm2        = AssociateSwarm(Swarm2,W,Z);
                    Cluster       = ClusterDivision(Swarm2,numClusters,clusterSize);
                end
                if stage == -1 && generation >= transitionGen+mergeCount*mergeInterval && mergeCount~=4
                    mergeCount  = mergeCount+1;
                    clusterSize = clusterSize*2;
                    [W,~]       = UniformPoint(floor(Problem.N/clusterSize),Problem.M);
                    numClusters = size(W,1);
                    Cluster     = ClusterDivision(Swarm2,numClusters,clusterSize);
                end
                if stage == 1
                    epsilon = inf;
                else
                    epsilon = 0;
                end
                Rank              = randperm(length(Swarm1),floor(length(Swarm1)/2)*2);
                Loser             = Rank(1:end/2);
                Winner            = Rank(end/2+1:end);
                Change            = Fitness1(Loser) <= Fitness1(Winner);
                Temp              = Winner(Change);
                Winner(Change)    = Loser(Change);
                Loser(Change)     = Temp;
                OffspringL        = CompetitiveLearning(Problem,Swarm1(Loser),Swarm1(Winner));
                Z                 = min([Z;OffspringL.objs],[],1);
                ParentW           = Swarm1(Winner);
                WinnerFitness     = CalFitness(ParentW.objs,ParentW.cons,epsilon);
                MatingPool        = TournamentSelection(2,Problem.N/2,WinnerFitness);
                OffspringW        = OperatorGA(Problem,ParentW(MatingPool));
                Z                 = min([Z;OffspringW.objs],[],1);
                NewCluster        = UnbalancedUpdate(Problem,Cluster,numClusters,clusterSize,W,Z);
                Offspring2        = [NewCluster{:}];
                [Swarm1,Fitness1] = SwarmEnvironmentalSelection([Swarm1 OffspringL OffspringW Offspring2],Problem.N,stage == -1);
                for i = 1 : numClusters
                    Cluster{i} = [Cluster{i} NewCluster{i}];
                end
                Cluster    = ClusterEnvironmentalSelection(Cluster,[OffspringL OffspringW],W,Z);
                Swarm2     = [Cluster{:}];
                generation = ceil(Problem.FE/(Problem.N*2));
                if ~mod(ceil(Problem.FE/(Problem.N*2)),ceil(0.1*Problem.maxFE/(Problem.N*2)))
                    Sw = [Swarm2,Swarm1];
                    W  = ReferenceVectorAdaptation(Sw.objs,W);
                end
                Z       = min([Z;Swarm2.objs],[],1);
                Z       = min([Z;Swarm1.objs],[],1);
                Archive = ArchiveUpdate([Archive Swarm2 Swarm1],Problem.N);
                if Problem.FE >= Problem.maxFE
                    Swarm2 = Archive;
                end
            end
        end
    end
end

function maxChange = MaximumChange(IdealHistory,NadirHistory,generation,historyLength)
    delta_value = 1e-6 * ones(1,size(IdealHistory,2));
    rz          = abs((IdealHistory(generation,:) - IdealHistory(generation - historyLength + 1,:)) ./ max(IdealHistory(generation - historyLength + 1,:),delta_value));
    nrz         = abs((NadirHistory(generation,:) - NadirHistory(generation - historyLength + 1,:)) ./ max(NadirHistory(generation - historyLength + 1,:),delta_value));
    maxChange   = max([rz, nrz]);
end

function NewSwarm = AssociateSwarm(Swarm,W,Z)
    N        = size(Swarm,2);
    NewSwarm = [];
    for i = 1 : N/2
        [~,Rank] = sort(acos(1-pdist2(W(i,:),(Swarm.objs-repmat(Z,N,1)),'cosine')));
        Next     = false(N,1);
        Next(Rank(1:2)) = true;
        NewSwarm = [NewSwarm Swarm(Next)];
    end
end

function Cluster = ClusterDivision(Swarm,numClusters,clusterSize)
    Cluster = cell(numClusters,1);
    N       = size(Swarm,2);
    if mod(N,numClusters) == 0
        for i = 1 : numClusters
            Cluster{i} = Swarm((i-1)*clusterSize+1:i*clusterSize);
        end
    else
        for i = 1 : numClusters
            if i ~= numClusters
                Cluster{i} = Swarm((i-1)*clusterSize+1:i*clusterSize);
            else
                Cluster{i} = Swarm((i-1)*clusterSize+1:N);
            end
        end
    end
end

function V = ReferenceVectorAdaptation(PopObj,V)
    A = max(PopObj,[],1)-min(PopObj,[],1);
    if ~all(A==0)
        V = V.*repmat(max(PopObj,[],1)-min(PopObj,[],1),size(V,1),1);
    end
end

function Offspring = CompetitiveLearning(Problem,Loser,Winner)
    LoserDec  = Loser.decs;
    WinnerDec = Winner.decs;
    [N,D]     = size(LoserDec);
    LoserVel  = Loser.adds(zeros(N,D));
    Winner.adds(zeros(N,D)); 
    r1        = repmat(rand(N,1),1,D);
    r2        = repmat(rand(N,1),1,D);
    OffVel    = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
    OffDec    = LoserDec + OffVel + r1.*(OffVel-LoserVel);
    OffDec    = PolynomialMutation(Problem,OffDec);
    Offspring = Problem.Evaluation(OffDec,OffVel);
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