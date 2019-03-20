function LMOCSO(Global)
% <algorithm> <L>
% Large-scale multi-objective competitive swarm optimization algorithm

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zheng, X. Zhang, and Y. Jin, Efficient large-scale
% multi-objective optimization based on a competitive swarm optimizer, IEEE
% Transactions on Cybernetics, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    [V,Global.N] = UniformPoint(Global.N,Global.M);
    Population   = Global.Initialization();
    Population   = EnvironmentalSelection(Population,V,(Global.gen/Global.maxgen)^2);
    
    %% Optimization
    while Global.NotTermination(Population)
        Fitness = calFitness(Population.objs);
        if length(Population) >= 2
            Rank = randperm(length(Population),floor(length(Population)/2)*2);
        else
            Rank = [1,1];
        end
        Loser  = Rank(1:end/2);
        Winner = Rank(end/2+1:end);
        Change = Fitness(Loser) >= Fitness(Winner);
        Temp   = Winner(Change);
        Winner(Change) = Loser(Change);
        Loser(Change)  = Temp;
        Offspring      = Operator(Population(Loser),Population(Winner));
        Population     = EnvironmentalSelection([Population,Offspring],V,(Global.gen/Global.maxgen)^2);
    end
end

function Fitness = calFitness(PopObj)
% Calculate the fitness by shift-based density

    N      = size(PopObj,1);
    fmax   = max(PopObj,[],1);
    fmin   = min(PopObj,[],1);
    PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    Dis    = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    Fitness = min(Dis,[],2);
end