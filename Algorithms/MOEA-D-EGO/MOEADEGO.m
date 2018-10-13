function MOEADEGO(Global)
% <algorithm> <H-N>
% Expensive Multiobjective Optimization by MOEA/D with Gaussian Process
% Model
% Ke ---  5 --- The number of function evaluations at each generation
% delta --- 0.9 --- The probability of choosing parents locally
% nr ---  2 --- Maximum number of solutions replaced by each offspring
% L1 --- 80 --- The maximal number of points used for building a local model
% L2 --- 20 --- The maximal number of points used for building a local model
% operator  --- DE

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    %% Parameter setting
    [Ke,delta,nr,L1,L2] = Global.ParameterSet(5,0.9,2,80,20);
    
	%% Generate random population
	NI = 11*Global.D-1;
    P  = lhsamp(NI,Global.D);
    Population = INDIVIDUAL(repmat(Global.upper-Global.lower,NI,1).*P+repmat(Global.lower,NI,1));

    %% Optimization
    while Global.NotTermination(Population)
        PopDec     = EGOSelect(Global,Population,L1,L2,Ke,delta,nr);
        Offspring  = INDIVIDUAL(PopDec);
        Population = [Population,Offspring];
    end
end