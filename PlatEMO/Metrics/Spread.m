function score = Spread(Population,optimum)
% <min> <multi/many> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <multimodal/none> <sparse/none> <dynamic/none>
% Spread

%------------------------------- Reference --------------------------------
% Y. Wang, L. Wu, and X. Yuan, Multi-objective self-adaptive differential
% evolution with elitist archive and crowding entropy-based diversity
% measure, Soft Computing, 2010, 14(3): 193-209.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.best.objs;
    if size(PopObj,2) ~= size(optimum,2)
        score = nan;
    else
        Dis1  = pdist2(PopObj,PopObj);
        Dis1(logical(eye(size(Dis1,1)))) = inf;
        [~,E] = max(optimum,[],1);
        Dis2  = pdist2(optimum(E,:),PopObj);
        d1    = sum(min(Dis2,[],2));
        d2    = mean(min(Dis1,[],2));
        score = (d1+sum(abs(min(Dis1,[],2)-d2))) / (d1+(size(PopObj,1)-size(PopObj,2))*d2);
    end
end