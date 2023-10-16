function P = BSPTreeLearning(P,best,lambda)
% Binary space partition tree based learning strategy

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    index = rand(1,size(P,1)) < lambda;
    P(index,1:best.level-1) = repmat(best.dec(1:best.level-1),sum(index),1);
end