function score = RVMM_IGD(PopObj,optimum)
% Inverted generational distance

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if size(PopObj,2) ~= size(optimum,2)
        score = nan;
    else
        score = mean(min(pdist2(optimum,PopObj),[],2));
    end
end