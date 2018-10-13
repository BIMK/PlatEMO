function Score = Coverage(PopObj,PF)
% <metric> <min>
% Coverage

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Domi = false(1,size(PopObj,1));
    for i = 1 : size(PF,1)
        Domi(sum(repmat(PF(i,:),size(PopObj,1),1)-PopObj<=0,2)==size(PopObj,2)) = true;
    end
    Score = sum(Domi) / size(PopObj,1);
end