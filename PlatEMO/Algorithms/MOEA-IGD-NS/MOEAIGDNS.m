function MOEAIGDNS(Global)
% <algorithm> <M>
% Multi-objective evolutionary algorithm based on an enhanced IGD

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zhang, R. Cheng, and Y. Jin, A multi-objective evolutionary
% algorithm based on an enhanced inverted generational distance metric,
% Proceedings of the IEEE Congress on Evolutionary Computation, 2016,
% 5222-5229.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the sampling points and random population
    Population = Global.Initialization();
    Archive    = UpdateArchive(Population,5*Global.N);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = randi(Global.N,1,Global.N);
        Offspring  = GA(Population(MatingPool));
        Archive    = UpdateArchive([Archive,Offspring],5*Global.N);
        Population = EnvironmentalSelection([Population,Offspring],Archive.objs,Global.N);
    end
end