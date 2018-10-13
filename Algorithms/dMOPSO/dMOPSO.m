function dMOPSO(Global)
% <algorithm> <A-G>
% A Multi-objective Particle Swarm Optimizer Based on Decomposition
% Ta --- 2 --- Age threshold
% operator --- dMOPSO_operator

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Ta = Global.ParameterSet(2);

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    
    %% Generate random population
    Population = Global.Initialization();
    Age        = zeros(1,Global.N);             % Age of each particle
    Z          = min(Population.objs,[],1);     % Ideal point
    Pbest      = Population;                    % Personal best
    Gbest      = UpdateGbest(W,Population,Z);   % Global best

    %% Optimization
    while Global.NotTermination(Gbest)
        % Generate offsprings
        Population(Age<Ta) = Global.Variation([Population(Age<Ta),Pbest(Age<Ta),Gbest(Age<Ta)],inf,@dMOPSO_operator);
        for i = find(Age>=Ta)
            temp          = randn(1,Global.D);
            Population(i) = INDIVIDUAL((Gbest(i).dec-Pbest(i).dec)/2+temp.*abs(Gbest(i).dec-Pbest(i).dec));
            Age(i)        = 0;
        end
        
        % Update Z
        Z = min(Z,min(Population.objs,[],1));
        
        % Update personal best
        PbeObj    = Pbest.objs - repmat(Z,Global.N,1);
        PopObj    = Population.objs - repmat(Z,Global.N,1);
        normW     = sqrt(sum(W.^2,2));
        normPbe   = sqrt(sum(PbeObj.^2,2));
        normPop   = sqrt(sum(PopObj.^2,2));
        CosinePbe = sum(PbeObj.*W,2)./normW./normPbe;
        CosinePop = sum(PopObj.*W,2)./normW./normPop;
        g_old     = normPbe.*CosinePbe + 5*normPbe.*sqrt(1-CosinePbe.^2);
        g_new     = normPop.*CosinePop + 5*normPop.*sqrt(1-CosinePop.^2);
        Pbest(g_new<=g_old) = Population(g_new<=g_old);
        Age(g_new<=g_old)   = 0;
        Age(g_new>g_old)    = Age(g_new>g_old) + 1;
        
        % Update global best
        Gbest = UpdateGbest(W,[Gbest,Population],Z);
    end
end