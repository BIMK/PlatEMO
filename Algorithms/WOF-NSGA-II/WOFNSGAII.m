function WOFNSGAII(Global)
% <algorithm> <O-Z>
% A Framework for Large-scale Multi-objective Optimization based on Problem
% Transformation
% t1    --- 1000 --- Number of evaluations for each optimization of original problem
% t2    ---  500 --- Number of evaluations for each optimization of transformed problem
% q     ---    3 --- Number of chosen solutions
% gamma ---    4 --- Number of groups
% delta ---  0.5 --- Ratio of evaluations for optimization of both original and transformed problems

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [t1,t2,q,gamma,delta] = Global.ParameterSet(1000,500,3,4,0.5);

    %% Generate random population
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);

    %% Optimization of both original and transformed problems
    while Global.evaluated <= delta*Global.evaluation
        % Optimization of original problem
        nEvaluated = Global.evaluated;
        while Global.NotTermination(Population) && Global.evaluated <= nEvaluated + t1
            MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
            Offspring  = Global.Variation(Population(MatingPool));
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
        end
        % Optimization of transformed problem
        [~,rank] = sortrows([FrontNo',-CrowdDis']);
        for k = 1 : q
            W = WeightingOptimization(Population(rank(k)).dec,t2,gamma,Global.lower,Global.upper);
            Population = [Population,INDIVIDUAL(Population.decs.*repmat(W,length(Population),1))];
        end
        [~,noDuplicate] = unique(Population.objs,'rows');
        Population      = Population(noDuplicate);
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,min(Global.N,length(Population)));
    end
    
    %% Optimization only original problem
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end