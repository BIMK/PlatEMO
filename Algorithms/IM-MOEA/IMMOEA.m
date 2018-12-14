function IMMOEA(Global)
% <algorithm> <I>
% Inverse modeling based multiobjective evolutionary algorithm
% K --- 10 --- Number of reference vectors

%------------------------------- Reference --------------------------------
% R. Cheng, Y. Jin, K. Narukawa, and B. Sendhoff, A multiobjective
% evolutionary algorithm using Gaussian process-based inverse modeling,
% IEEE Transactions on Evolutionary Computation, 2015, 19(6): 838-856.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    K = Global.ParameterSet(10);

    %% Generate random population
    [W,K] = UniformPoint(K,Global.M);
    W     = fliplr(sortrows(fliplr(W)));
    Global.N      = ceil(Global.N/K)*K;
    Population    = Global.Initialization();
    [~,partition] = max(1-pdist2(Population.objs,W,'cosine'),[],2);
    
    %% Optimization
    while Global.NotTermination(Population)
        % Modeling and reproduction
        for k = unique(partition)'
            Population = [Population,Operator(Population(partition==k))];
        end
        % Environmental selection
        [~,partition] = max(1-pdist2(Population.objs,W,'cosine'),[],2);
        for k = unique(partition)'
            current = find(partition==k);
            if length(current) > Global.N/K
                Del = EnvironmentalSelection(Population(current),Global.N/K);
                Population(current(Del)) = [];
                partition(current(Del))  = [];
            end
        end
    end
end