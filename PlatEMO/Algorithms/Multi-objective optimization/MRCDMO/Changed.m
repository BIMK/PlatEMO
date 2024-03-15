function changed = Changed(Problem,Population)
% Detect whether the problem changes

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    RePop1  = Population(randperm(end,ceil(end/10)));
    RePop2  = Problem.Evaluation(RePop1.decs);
    changed = ~isequal(RePop1.objs,RePop2.objs) || ~isequal(RePop1.cons,RePop2.cons);
end