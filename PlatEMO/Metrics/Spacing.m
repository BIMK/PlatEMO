function score = Spacing(Population,~)
% <min> <multi/many> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <multimodal/none> <sparse/none> <dynamic/none>
% Spacing

%------------------------------- Reference --------------------------------
% J. R. Schott, Fault tolerant design using single and multicriteria
% genetic algorithm optimization, Master's thesis, Department of
% Aeronautics and Astronautics, Massachusetts Institute of Technology,
% 1995.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.best.objs;
    if isempty(PopObj)
        score = nan;
    else
        Distance = pdist2(PopObj,PopObj,'cityblock');
        Distance(logical(eye(size(Distance,1)))) = inf;
        score = std(min(Distance,[],2));
    end
end