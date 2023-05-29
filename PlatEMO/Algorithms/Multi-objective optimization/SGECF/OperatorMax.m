function [OffDec,OffMask] = OperatorMax(Problem,winpop,windec,winmask,winrank,SubPopulation,SubDec,SubMask,Rank,Fitness,thetamid,REAL)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    winPool     = TournamentSelection(2,size(SubPopulation,2),winrank);
    Parent1Dec  = windec(winPool,:);
    Parent1Mask = winmask(winPool,:);

    MatingPool  = TournamentSelection(2,size(SubPopulation,2),Rank);
    Parent2Dec  = SubDec(MatingPool,:);
    Parent2Mask = SubMask(MatingPool,:);

    ParentDec = [Parent1Dec;Parent2Dec];
    [N,D]     = size(ParentDec);

    %% Perform uniform crossover and bit-flip mutation for Mask
    % Uniform crossover
    k = rand(N/2,D) < 0.5;
    k(repmat(rand(N/2,1)>1,1,D)) = false;
    OffMask    = Parent1Mask;
    OffMask(k) = Parent2Mask(k);

    %% Mutation for mask
    for i = 1 : N/2
        if sum(OffMask(i,:))<=thetamid
            % Bit-flip mutation
            Site = rand(1,D) < 1/D;
            OffMask(i,Site) = ~OffMask(i,Site);
        else
            index1 = find(OffMask(i,:));
            index0 = find(~OffMask(i,:));
            maxN = floor(size(index1,2) - thetamid);
            index = index1(TSmin(-Fitness(index1),maxN));
            OffMask(i,index) = 0;
        end
    end

    %% Crossover and mutation for dec
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