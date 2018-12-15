function Score = Spread(PopObj,PF)
% <metric> <min>
% Spread

%------------------------------- Reference --------------------------------
% Y. Wang, L. Wu, and X. Yuan, Multi-objective self-adaptive differential
% evolution with elitist archive and crowding entropy-based diversity
% measure, Soft Computing, 2010, 14(3): 193-209.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Dis1  = pdist2(PopObj,PopObj);
    Dis1(logical(eye(size(Dis1,1)))) = inf;
    [~,E] = max(PF,[],1);Draw(PF(E,:),'*b');PF(E,:)
    Dis2  = pdist2(PF(E,:),PopObj);
    d1    = sum(min(Dis2,[],2));
    d2    = mean(min(Dis1,[],2));
    Score = (d1+sum(abs(min(Dis1,[],2)-d2))) / (d1+(size(PopObj,1)-size(PopObj,2))*d2);
end