function State = EstimateState(Archive,StateNum)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    State    = zeros(1, StateNum);
    Cons     = sum(max(Archive.cons,0),2);
    FIndex   = Cons == 0;
    State(1) = mean(Cons);
    if any(FIndex)
        [DNum,NNum] = Detect(Archive, FIndex);
        State(2)    = DNum/numel(Archive);
        State(3)    = NNum/numel(Archive);
    end    
end

function [DNum, NNum] = Detect(Archive, FIndex)
    DNum = 0;
    NNum = 0;
    FSet = Archive(FIndex);
    FNum = numel(FSet);
    ISet = Archive(~FIndex);
    INum = numel(ISet);
    for i = 1 : INum
        Del       = all(repmat(ISet(i).objs, FNum, 1) < FSet.objs, 2);
        FSet(Del) = [];
        DNum      = DNum + sum(Del);
        FNum      = numel(FSet);
        if FNum == 0
            return;
        end
    end
    NNum = sum(FIndex) - DNum;
end