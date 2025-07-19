function [Population,FrontNo,CrowdDis] = Reinitialization(Problem,Population,zeta)
% Re-initialize solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N          = floor(length(Population)*zeta/2)*2;
    Selected   = randperm(length(Population),N);
    Population(Selected)   = Problem.Initialization(N);
    unSelected = setdiff(1:length(Population),Selected);
    Population(unSelected) = Problem.Evaluation(Population(unSelected).decs);
    [~,FrontNo,CrowdDis]   = EnvironmentalSelection(Population,length(Population));
end