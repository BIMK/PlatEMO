function eMOEA(Global)
% <algorithm> <E>
% ¦Å multi-objective evolutionary algorithm
% epsilon --- 0.06 --- The parameter in grid location calculation

%------------------------------- Reference --------------------------------
% K. Deb, M. Mohan, and S. Mishra, Towards a quick computation of
% well-spread Pareto-optimal solutions, Proceedings of the International
% Conference on Evolutionary Multi-Criterion Optimization, 2003, 222-236.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
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
            Offspring  = GAhalf([Population(p),Archive(q)]);
            Population = UpdatePopulation(Population,Offspring);
            Archive    = UpdateArchive(Archive,Offspring,epsilon);
        end
    end
end