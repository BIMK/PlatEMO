function [better_Cpop,bad_Cpop] = DimSelect(Population,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yingwei Li

    PopObj          = Population.objs;
    [Convergence,~] = CalFitness(PopObj);
    [~,Cidx]        = sort(Convergence);
    better_Cpop     = Population(Convergence==0);
    temple          = length(Convergence)-N+1;
    bad_Cpop        = Population(Cidx(temple:end));
end