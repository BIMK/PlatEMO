function MaOEAIGD(Global)
% <algorithm> <M>
% IGD based many-objective evolutionary algorithm
% DNPE --- --- Number of evaluations for nadir point estimation

%------------------------------- Reference --------------------------------
% Y. Sun, G. G. Yen, and Z. Yi, IGD indicator-based evolutionary algorithm
% for many-objective optimization problems, IEEE Transactions on
% Evolutionary Computation, 2018.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    DNPE = Global.ParameterSet(100*Global.N);

    %% Nadir point estimation
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    % Optimize Global.M single-objective optimization problems
    Population   = Global.Initialization();
    while Global.NotTermination(Population) && Global.evaluated < DNPE
        Offspring  = GA(Population(randi(end,1,Global.N)),{0.9,20,1,20});
        Population = [Population,Offspring];
        [~,rank]   = sort(Fitness(Population.objs),1);
        Population = Population(unique(rank(1:ceil(Global.N/Global.M),:)));
    end
    % Find the nadir point and ideal point
    [~,ext] = min(Fitness(Population.objs),[],1);
    zmax    = diag(Population(ext).objs)';
    zmin    = min(Population.objs,[],1);
    zmax(zmax<1e-6) = 1;
    % Transform the points into the Utopian Pareto front
    W = W.*repmat(zmax-zmin,Global.N,1) + repmat(zmin,Global.N,1);

    %% Generate random population
    Population = Global.Initialization();
    [Population,Rank,Dis] = EnvironmentalSelection(Population,W,Global.N);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,Rank,min(Dis,[],2));
        Offspring  = GAhalf(Population(MatingPool));
        [Population,Rank,Dis] = EnvironmentalSelection([Population,Offspring],W,Global.N);
    end
end

function fit = Fitness(PopObj)
% Calculate the objective value of each solution on each single-objective
% optimization problem in nadir point estimation

    fit   = zeros(size(PopObj));
    for i = 1 : size(PopObj,2)
        fit(:,i) = abs(PopObj(:,i)) + 100*sum(PopObj(:,[1:i-1,i+1:end]).^2,2);
    end
end