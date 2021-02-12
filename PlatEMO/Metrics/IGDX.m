function score = IGDX(Population,optimum)
% <min>
% Inverted generational distance in the decision space

%------------------------------- Reference --------------------------------
% A. Zhou, Q. Zhang, and Y. Jin, Approximating the set of Pareto-optimal
% solutions in both the decision and objective spaces by an estimation of
% distribution algorithm, IEEE Transactions on Evolutionary Computation,
% 2009, 13(5): 1167-1189.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopDec = Population.decs;
    if size(PopDec,2) ~= size(optimum,2)
        score = nan;
    else
        score = mean(min(pdist2(optimum,PopDec),[],2));
    end
end