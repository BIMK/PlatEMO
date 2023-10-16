function Fitness = CalFitness(C,Population)
% Calculate the fitness of each solution in terms of a single level

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.objs;
    PopCon = Population.cons;
    if any(isnan(PopObj(:,1)))  % Lower level
        PopObj = PopObj(:,2);
        PopCon = PopCon(:,C+1:end);
    else                        % Upper level
        PopObj = PopObj(:,1);
        PopCon = PopCon(:,1:C);
    end
    if isempty(PopCon)
        PopCon = zeros(size(PopObj,1),1);
    else
        PopCon = sum(max(0,PopCon),2);
    end
    Feasible = PopCon <= 0;
    Fitness  = Feasible.*PopObj + ~Feasible.*(PopCon+1e10);
end