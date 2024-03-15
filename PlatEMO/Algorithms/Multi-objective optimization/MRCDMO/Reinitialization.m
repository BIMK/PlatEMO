function [Population,FrontNo,CrowdDis] = Reinitialization(Problem,Population,type,zeta)
% Re-initialize solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = floor(length(Population)*zeta/2)*2;
    Selected = randperm(length(Population),N);
    if type == 1
        Population(Selected) = OperatorGA(Problem,Population(Selected),{0,0,inf,10});
    else
        Population(Selected) = Problem.Initialization(N);
    end
    unSelected = setdiff(1:length(Population),Selected);
    Population(unSelected) = Problem.Evaluation(Population(unSelected).decs);
    [~,FrontNo,CrowdDis]   = EnvironmentalSelection(Population,length(Population));
end