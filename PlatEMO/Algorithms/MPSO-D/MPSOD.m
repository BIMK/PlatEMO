function MPSOD(Global)
% <algorithm> <M>
% Multi-objective particle swarm optimization algorithm based on
% decomposition

%------------------------------- Reference --------------------------------
% C. Dai, Y. Wang, and M. Ye, A new multi-objective particle swarm
% optimization algorithm based on decomposition, Information Sciences,
% 2015, 325: 541-557.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
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
        [Parent,Pbest,Gbest] = MatingSelection(Population.objs,B,W,Z);
        Offspring  = Operator(Population(Parent),Population(Pbest),Population(Gbest));
        Z          = min([Z;Offspring.objs],[],1);
        Population = Classification(Global,[Population,Offspring],W,Z);
    end
end