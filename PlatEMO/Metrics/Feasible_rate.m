function score = Feasible_rate(Population,~)
% <max> <single/multi/many> <real/integer/label/binary/permutation> <large/none> <constrained> <expensive/none> <multimodal/none> <sparse/none> <dynamic/none> <multitask/none> <bilevel/none> <robust/none>
% The rate of feasible solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	score = mean(all(Population.cons<=0,2));
end