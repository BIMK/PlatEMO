function AGEMOEA(Global)
% <algorithm> <H-N>
% An Adaptive Evolutionary Algorithm based on Non-Euclidean Geometry for Many-objective Optimization
% submitted to GECCO 2019

    %% Generate random population
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    
    %% Optimization
    while Global.NotTermination(Population)
      MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
      Offspring  = GA(Population(MatingPool));
      [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end