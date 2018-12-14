function SRA(Global)
% <algorithm> <S>
% Stochastic ranking algorithm

%------------------------------- Reference --------------------------------
% B. Li, K. Tang, J. Li, and X. Yao, Stochastic ranking algorithm for
% many-objective optimization based on multiple indicators, IEEE
% Transactions on Evolutionary Computation, 2016, 20(6): 924-938.
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
        MatingPool = randi(Global.N,1,2*Global.N);
        Offspring  = GAhalf(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],Global.N,unifrnd(0.4,0.6));
    end
end