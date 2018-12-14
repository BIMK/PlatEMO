function CMOPSO(Global) 
% <algorithm> <C>
% Competitive mechanism based multi-objective particle swarm optimizer

%------------------------------- Reference --------------------------------
% X. Zhang, X. Zheng, R. Cheng, J. Qiu, and Y. Jin, A competitive mechanism
% based multi-objective particle swarm optimizer with fast convergence,
% Information Sciences, 2018, 427: 63-76.
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
        Offspring  = Operator(Population);
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end