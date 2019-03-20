function LSMOF(Global)
% <algorithm> <L>
% Large-scale multi-objective optimization framework with NSGA-II
% G1       --- 10 --- The generation of weight optimization with DE
% SubN     --- 30 --- The population size of the transferred problem
% operator ---  1 --- Original operators 1. GA 2. DE

%------------------------------- Reference --------------------------------
% C. He, L. Li, Y. Tian, X. Zhang, R. Cheng, Y. Jin, and X. Yao,
% Accelerating large-scale multi-objective optimization via problem
% reformulation, IEEE Transactions on Evolutionary Computation, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    %% Parameter settings
    [wD,SubN,Operator] = Global.ParameterSet(10,30,2);
    Population         = Global.Initialization();
    G = floor(Global.evaluation*0.05/(SubN*2*wD));
    
    %% Optimization
    while Global.NotTermination(Population)
        if Global.evaluated < 0.6*Global.evaluation   
            Archive    = WeightOptimization(Global,G,Population,wD,SubN);
            Population = EnvironmentalSelection([Population,Archive],Global.N);
        else
            Population = subNSGAII(Population,Operator,Global.N);
        end
    end
end