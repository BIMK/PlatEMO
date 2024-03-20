function score = Upper_level_Min_value(Population,~)
% <min> <multi> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <sparse/none> <bilevel>
% The minimum objective value of the upper level (for bilevel optimization)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Feasible = find(all(Population.cons<=0,2));
    if ~isempty(Feasible)
        PopObj = Population(Feasible).objs;
        score  = min(PopObj(:,1));
    else
        score = nan;
    end
end