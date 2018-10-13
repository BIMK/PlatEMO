function eMOEA(Global)
% <algorithm> <A-G>
% Towards a Quick Computation of Well-Spread Pareto-Optimal Solutions
% epsilon --- 0.06 --- The parameter in grid location calculation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    epsilon = Global.ParameterSet(0.06);

    %% Generate random population
    Population = Global.Initialization();
    PopGrid    = floor((Population.objs-repmat(min(Population.objs,[],1),Global.N,1))/epsilon);
    eFrontNO   = NDSort(PopGrid,1);
    Archive    = Population(eFrontNO==1);
    
    %% Optimization
    while Global.NotTermination(Archive)
        for i = 1 : Global.N
            k    = randperm(Global.N,2);
            Domi = any(Population(k(1)).obj<Population(k(2)).obj) - any(Population(k(1)).obj>Population(k(2)).obj);
            p    = k((Domi==-1)+1);
            q    = randi(length(Archive));
            Offspring  = Global.Variation([Population(p),Archive(q)],1);
            Population = UpdatePopulation(Population,Offspring);
            Archive    = UpdateArchive(Archive,Offspring,epsilon);
        end
    end
end