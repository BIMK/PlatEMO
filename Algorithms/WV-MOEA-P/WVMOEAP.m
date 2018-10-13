function WVMOEAP(Global)
% <algorithm> <O-Z>
% A Weight Vector Based Multi-Objective Optimization Algorithm with
% Preference (Chinese)
% Points ---      --- Set of preferred points
% b      --- 0.05 --- Extent of preferred region
% operator        --- DE

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [Points,b] = Global.ParameterSet(ones(1,Global.M),0.05);

    %% Generate the weight vectors and random population
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T = ceil(Global.N/10);

    %% Map the weight vectors
    Dis   = pdist2(W,W);
    B     = zeros(Global.N,T);
    Group = ceil((1:Global.N)/Global.N*size(Points,1));
    for i = unique(Group)
        % The weight vectors around the i-th preferred point
        Current = find(Group==i);
        % Map the weight vectors
        W(Current,:) = 2*b.*W(Current,:) + Points(i,:) + b;
        % Detect the neighbours of each vector
        [~,rank]     = sort(Dis(Current,Current),2);
        B(Current,:) = Current(rank(:,1:T));
    end

    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        % For each group
        for i = unique(Group)
            Current = find(Group==i);
            % Generate an offspring for each solution
            P = zeros(length(Current),2);
            for j = 1 : size(P,1)
                if rand < 0.9
                    P(j,:) = B(j,randperm(size(B,2),2));
                else
                    P(j,:) = Current(randperm(length(Current),2));
                end
            end
            Offspring = Global.Variation(Population([Current,P(:)']),inf,@DE);
            % Environmental selection
            Z = min(Z,min(Offspring.objs,[],1));
            Population(Current) = STM([Population(Current),Offspring],W(Current,:),Z,Z+ones(1,Global.M));
        end
    end
end