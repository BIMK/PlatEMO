function score = Min_value(Population,~)
% <min> <single> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <sparse/none> <dynamic/none>
% The minimum objective value

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	score = Population.best.objs;
    if isempty(score)
        score = nan;
    end
end