function MPSOD(Global)
% <algorithm> <H-N>
% A new multi-objective particle swarm optimization algorithm based on
% decomposition
% operator --- MPSOD_operator

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    W = W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));
    T = ceil(Global.N/10);
    
	%% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %% Generate random population
    Population = Global.Initialization(2*Global.N);
    Z          = min(Population.objs,[],1);
    Population = Classification(Global,Population,W,Z);
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,B,W,Z);
        Offspring  = Global.Variation(Population(MatingPool),inf,@MPSOD_operator);
        Z          = min([Z;Offspring.objs],[],1);
        Population = Classification(Global,[Population,Offspring],W,Z);
    end
end