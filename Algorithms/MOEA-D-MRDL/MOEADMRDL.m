function MOEADMRDL(Global)
% <algorithm> <H-N>
% Online Diversity Assessment in Evolutionary Multiobjective Optimization:
% A Geometrical Perspective
% gamma --- 20 --- Maximum allowable maximum relative diversity loss
% nmov  --- 10 --- Size of moving average
% operator     --- MOEADMRDL_operator

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [gamma,nmov] = Global.ParameterSet(20,10);
    
    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);     % Weight vectors
    T = ceil(Global.N/10);                              % Neighbourhood size
    
    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T)';                          % The neighbours of each weight vector
    
    %% Generate random population
    Population = Global.Initialization();   % The initial population
    Z = min(Population.objs,[],1);          % The ideal point

    %% Optimization
    disC      = 20;     % The distribution index of SBX
    Pn        = 0;      % The standard deviation of Gaussian permutation
    allEgamma = [];     % The MRDL array storing all the past MRDL values
    while Global.NotTermination(Population)
        % Select parents
        MatingPool = [B(randi(T,1,Global.N)+(0:Global.N-1)*T),B(randi(T,1,Global.N)+(0:Global.N-1)*T)];
        % Generate offsprings
        Offspring  = Global.Variation(Population(MatingPool),Global.N,@MOEADMRDL_operator,{1,disC,1,20,Pn});
        % Update the ideal point
        Z          = min([Z;Offspring.objs],[],1);
        % Environmental selection
        [Population,Egamma] = EnvironmentalSelection(Population,Offspring,W,B',Z,gamma);
        % Operator adaption
        [disC,Pn,allEgamma] = OperatorAdaption(disC,Pn,allEgamma,Egamma,nmov);
    end
end