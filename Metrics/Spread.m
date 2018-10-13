function Score = Spread(PopObj,PF)
% <metric> <min>
% Spread

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Dis1  = pdist2(PopObj,PopObj);
    Dis1(logical(eye(size(Dis1,1)))) = inf;
    [~,E] = max(PF,[],1);
    Dis2  = pdist2(PF(E,:),PopObj);
    d1    = sum(min(Dis2,[],2));
    d2    = mean(min(Dis1,[],2));
    Score = (d1+sum(abs(min(Dis1,[],2)-d2))) / (d1+(size(PopObj,1)-size(PopObj,2))*d2);
end