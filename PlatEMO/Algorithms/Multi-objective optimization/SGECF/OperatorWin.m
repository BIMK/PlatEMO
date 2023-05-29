function [OffDec,OffMask] = OperatorWin(Problem,SubPopulation,SubDec,SubMask,Rank,Fitness,REAL)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    MatingPool = TournamentSelection(2,2*size(SubPopulation,2),Rank);
    ParentDec = SubDec(MatingPool,:);
    ParentMask = SubMask(MatingPool,:);
    %% Parameter setting
    [N,D]       = size(ParentDec);
    randN = randperm(N);
    Parent1Dec  = ParentDec(randN(1:N/2),:);
    Parent2Dec  = ParentDec(randN(N/2+1:end),:);
    Parent1Mask = ParentMask(randN(1:N/2),:);
    Parent2Mask = ParentMask(randN(N/2+1:end),:);

    %% Perform uniform crossover and bit-flip mutation for Mask
    % Uniform crossover
    k = rand(N/2,D) < 0.5;
    k(repmat(rand(N/2,1)>1,1,D)) = false;
    OffMask    = Parent1Mask;
    OffMask(k) = Parent2Mask(k);
    % Bit-flip mutation
    Site = rand(N/2,D) < 1/D;
    OffMask(Site) = ~OffMask(Site);

    %% Crossover and mutation for dec
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