function score = DeltaP(Population,optimum)
% <min> <multi/many> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <multimodal/none> <sparse/none> <dynamic/none>
% Averaged Hausdorff distance

%------------------------------- Reference --------------------------------
% O. Schutze, X. Esquivel, A. Lara, and C. A. Coello Coello, Using the
% averaged Hausdorff distance as a performance measure in evolutionary
% multiobjective optimization, IEEE Transactions on Evolutionary
% Computation, 2012, 16(4): 504-522.
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
        Distance = pdist2(optimum,PopObj);
        score    = max(mean(min(Distance,[],1)),mean(min(Distance,[],2)));
    end
end