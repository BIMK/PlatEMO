function MatingPool = MatingStrategy(IMatrix)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiawei Yuan

Ps                   = 0.7;     % probability of partner selection
[Num,~]              = size(IMatrix);
[~,Neighboors]       = min(IMatrix,[],2);

AllInd               = 1:Num;
SpouseID             = Neighboors;
ChnageInd            = find(rand(1,Num)>Ps);
SpouseID(ChnageInd)  = AllInd(ceil(rand(1,length(ChnageInd))*Num));

AllInd               = reshape(AllInd,1,Num);
SpouseID             = reshape(SpouseID,1,Num);
MatingPool           = [AllInd,SpouseID];
end