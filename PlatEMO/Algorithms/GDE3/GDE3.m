function GDE3(Global)
% <algorithm> <G>
% Generalized differential evolution 3

%------------------------------- Reference --------------------------------
% S. Kukkonen and J. Lampinen, GDE3: The third evolution step of
% generalized differential evolution, Proceedings of the IEEE Congress on
% Evolutionary Computation, 2005, 443-450.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        Offspring  = DE(Population,Population(randi(Global.N,1,Global.N)),Population(randi(Global.N,1,Global.N)));
        Population = EnvironmentalSelection(Population,Offspring,Global.N);
    end
end