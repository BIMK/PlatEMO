function IMMOEA(Global)
% <algorithm> <H-N>
% A Multiobjective Evolutionary Algorithm using Gaussian Process based
% Inverse Modeling
% K --- 10 --- Number of reference vectors
% operator --- IMMOEA_operator

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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
            Population = [Population,Global.Variation(Population(partition==k),inf,@IMMOEA_operator)];
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