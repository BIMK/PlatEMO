function [OffDec,OffMask] = Operator(Population,Dec,Mask,Fitness,winIndex,loseIndex1,loseIndex2,Problem,thetamid,REAL)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    SubPopulation    = {};
    SubPopulation{1} = Population(winIndex);
    SubPopulation{2} = Population(loseIndex1);
    SubPopulation{3} = Population(loseIndex2);

    SubDec    = {};
    SubDec{1} = Dec(winIndex,:);
    SubDec{2} = Dec(loseIndex1,:);
    SubDec{3} = Dec(loseIndex2,:);

    SubMask    = {};
    SubMask{1} = Mask(winIndex,:);
    SubMask{2} = Mask(loseIndex1,:);
    SubMask{3} = Mask(loseIndex2,:);

    [SubPopulation,SubDec,SubMask,Rank] = InitRank(SubPopulation,SubDec,SubMask,Problem);
    [OffDec1,OffMask1] = OperatorWin(Problem,SubPopulation{1},SubDec{1},SubMask{1},Rank{1},Fitness,REAL);
    if size(SubPopulation{2},2)>0
        [OffDec2,OffMask2] = OperatorMin(Problem,SubPopulation{1},SubDec{1},SubMask{1},Rank{1},SubPopulation{2},SubDec{2},SubMask{2},Rank{2},Fitness,thetamid,REAL);
    else
        OffDec2  = [];
        OffMask2 = [];
    end
    if size(SubPopulation{3},2)>0
        [OffDec3,OffMask3] = OperatorMax(Problem,SubPopulation{1},SubDec{1},SubMask{1},Rank{1},SubPopulation{3},SubDec{3},SubMask{3},Rank{3},Fitness,thetamid,REAL);
    else
        OffDec3  = [];
        OffMask3 = [];
    end
    OffDec  = [OffDec1;OffDec2;OffDec3];
    OffMask = [OffMask1;OffMask2;OffMask3];
end