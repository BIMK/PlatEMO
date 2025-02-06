function [OffDec,OffMask] = Operator(Population,Dec,Mask,Fitness,winIndex,loseIndex1,loseIndex2,Problem,thetamid,REAL)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Divide population into three subpopulations
    SubPopulation = {};
    SubPopulation{1} = Population(winIndex);
    SubPopulation{2} = Population(loseIndex1);
    SubPopulation{3} = Population(loseIndex2);

    SubDec = {};
    SubDec{1} = Dec(winIndex,:);
    SubDec{2} = Dec(loseIndex1,:);
    SubDec{3} = Dec(loseIndex2,:);

    SubMask = {};
    SubMask{1} = Mask(winIndex,:);
    SubMask{2} = Mask(loseIndex1,:);
    SubMask{3} = Mask(loseIndex2,:);

    % Sort the solutions in population
    [SubPopulation,SubDec,SubMask,Rank] = InitRank(SubPopulation,SubDec,SubMask); 

    % Generate Offspring along the direction increasing sparsity to cosparsity
    if size(SubPopulation{2},2)>size(SubPopulation{3},2)&&size(SubPopulation{2},2)>floor(Problem.N/3)
        [OffDec1,OffMask1] = OperatorWin(Problem,[SubPopulation{1},SubPopulation{3}],[SubDec{1};SubDec{3}],[SubMask{1};SubMask{3}],[Rank{1},Rank{3}+max(Rank{1})],Fitness,REAL);
        [OffDec2,OffMask2] = OperatorMin(Problem,SubPopulation{1},SubDec{1},SubMask{1},Rank{1},SubPopulation{2},SubDec{2},SubMask{2},Rank{2},Fitness,thetamid,REAL);
        OffDec3 = [];
        OffMask3 = [];
    % Generate Offspring along the direction decreasing sparsity to cosparsity
    elseif size(SubPopulation{2},2)<size(SubPopulation{3},2)&&size(SubPopulation{3},2)>floor(Problem.N/3)
        [OffDec1,OffMask1] = OperatorWin(Problem,[SubPopulation{1},SubPopulation{2}],[SubDec{1};SubDec{2}],[SubMask{1};SubMask{2}],[Rank{1},Rank{2}+max(Rank{1})],Fitness,REAL);
        [OffDec3,OffMask3] = OperatorMax(Problem,SubPopulation{1},SubDec{1},SubMask{1},Rank{1},SubPopulation{3},SubDec{3},SubMask{3},Rank{3},Fitness,thetamid,REAL);
        OffDec2 = [];
        OffMask2 = [];
    else
        [OffDec1,OffMask1] = OperatorWin(Problem,SubPopulation{1},SubDec{1},SubMask{1},Rank{1},Fitness,REAL);
        if size(SubPopulation{2},2)>0
            [OffDec2,OffMask2] = OperatorMin(Problem,SubPopulation{1},SubDec{1},SubMask{1},Rank{1},SubPopulation{2},SubDec{2},SubMask{2},Rank{2},Fitness,thetamid,REAL);
        else
            OffDec2 = [];
            OffMask2 = [];
        end
        if size(SubPopulation{3},2)>0
            [OffDec3,OffMask3] = OperatorMax(Problem,SubPopulation{1},SubDec{1},SubMask{1},Rank{1},SubPopulation{3},SubDec{3},SubMask{3},Rank{3},Fitness,thetamid,REAL);
        else
            OffDec3 = [];
            OffMask3 = [];
        end
    end
    OffDec  = [OffDec1;OffDec2;OffDec3];
    OffMask = [OffMask1;OffMask2;OffMask3];
end

function [OffDec,OffMask] = OperatorWin(Problem,SubPopulation,SubDec,SubMask,Rank,Fitness,REAL)
    MatingPool = TournamentSelection(2,2*size(SubPopulation,2),Rank);
    ParentDec = SubDec(MatingPool,:);
    ParentMask = SubMask(MatingPool,:);
    % Parameter setting
    [N,D]       = size(ParentDec);
    randN = randperm(N);
    Parent1Dec  = ParentDec(randN(1:N/2),:);
    Parent2Dec  = ParentDec(randN(N/2+1:end),:);
    Parent1Mask = ParentMask(randN(1:N/2),:);
    Parent2Mask = ParentMask(randN(N/2+1:end),:);
    
    % Crossover
    k = rand(N/2,D) < 0.5;
    k(repmat(rand(N/2,1)>1,1,D)) = false;
    OffMask    = Parent1Mask;
    OffMask(k) = Parent2Mask(k);
    % Mutation
    Site = rand(N/2,D) < 1/D;
    OffMask(Site) = ~OffMask(Site);

    % Crossover and mutation for dec
    if REAL
        OffDec = OperatorGAhalf(Problem,ParentDec);
    else
        OffDec = ones(N/2,D);
    end
end

function index = TS(Fitness)
% Binary tournament selection

    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(2,1,Fitness);
    end
end

function [OffDec,OffMask] = OperatorMin(Problem,winpop,windec,winmask,winrank,SubPopulation,SubDec,SubMask,Rank,Fitness,thetamid,REAL)
    winPool     = TournamentSelection(2,size(SubPopulation,2),winrank);
    Parent1Dec  = windec(winPool,:);
    Parent1Mask = winmask(winPool,:);

    MatingPool  = TournamentSelection(2,size(SubPopulation,2),Rank);
    Parent2Dec  = SubDec(MatingPool,:);
    Parent2Mask = SubMask(MatingPool,:);

    ParentDec = [Parent1Dec;Parent2Dec];
    [N,D]     = size(ParentDec);
    
    % Crossover
    k = rand(N/2,D) < 0.5;
    k(repmat(rand(N/2,1)>1,1,D)) = false;
    OffMask    = Parent1Mask;
    OffMask(k) = Parent2Mask(k);

    % Mutation
    for i = 1 : N/2
        if sum(OffMask(i,:))>=thetamid
            % Bit-flip mutation
            Site = rand(1,D) < 1/D;
            OffMask(i,Site) = ~OffMask(i,Site);
        else
            index1 = find(OffMask(i,:));
            index0 = find(~OffMask(i,:));
            minN = floor(thetamid - size(index1,2));
            index = index0(TSmin(Fitness(index0),minN));
            OffMask(i,index) = 1;
        end
    end

    % Crossover and mutation for dec
    if REAL
        OffDec = OperatorGAhalf(Problem,ParentDec);
    else
        OffDec = ones(N/2,D);
    end
end

function index = TSmin(Fitness,minN)
% Binary tournament selection

    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(2,minN,Fitness);
    end
end

function [OffDec,OffMask] = OperatorMax(Problem,winpop,windec,winmask,winrank,SubPopulation,SubDec,SubMask,Rank,Fitness,thetamid,REAL)
    winPool     = TournamentSelection(2,size(SubPopulation,2),winrank);
    Parent1Dec  = windec(winPool,:);
    Parent1Mask = winmask(winPool,:);

    MatingPool  = TournamentSelection(2,size(SubPopulation,2),Rank);
    Parent2Dec  = SubDec(MatingPool,:);
    Parent2Mask = SubMask(MatingPool,:);

    ParentDec = [Parent1Dec;Parent2Dec];
    [N,D]     = size(ParentDec);
    
    % Crossover
    k = rand(N/2,D) < 0.5;
    k(repmat(rand(N/2,1)>1,1,D)) = false;
    OffMask    = Parent1Mask;
    OffMask(k) = Parent2Mask(k);

    % Mutation
    for i = 1 : N/2
        if sum(OffMask(i,:))<=thetamid
            % Bit-flip mutation
            Site = rand(1,D) < 1/D;
            OffMask(i,Site) = ~OffMask(i,Site);
        else
            index1 = find(OffMask(i,:));
            index0 = find(~OffMask(i,:));
            maxN = floor(size(index1,2) - thetamid);
            index = index1(TSmax(-Fitness(index1),maxN));
            OffMask(i,index) = 0;
        end
    end

    % Crossover and mutation for dec
    if REAL
        OffDec = OperatorGAhalf(Problem,ParentDec);
    else
        OffDec = ones(N/2,D);
    end
end

function index = TSmax(Fitness,minN)
% Binary tournament selection

    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(2,minN,Fitness);
    end
end