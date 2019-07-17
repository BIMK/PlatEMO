function AGEMOEA(Global)
% <algorithm> <A>
% Adaptive geometry estimation-based many-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% A. Panichella, An adaptive evolutionary algorithm based on non-euclidean
% geometry for many-objective optimization, Proceedings of the Genetic and
% Evolutionary Computation Conference, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Annibale Panichella

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