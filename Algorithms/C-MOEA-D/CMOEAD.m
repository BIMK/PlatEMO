function CMOEAD(Global)
% <algorithm> <A-G>
% An Evolutionary Many-Objective Optimization Algorithm Using
% Reference-Point Based Nondominated Sorting Approach, Part II: Handling
% Constraints and Extending to an Adaptive Approach

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
    T  = ceil(Global.N/10);
    nr = ceil(Global.N/100);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        % For each solution
        for i = 1 : Global.N      
            % Choose the parents
            if rand < 0.9
                P = B(i,randperm(size(B,2)));
            else
                P = randperm(Global.N);
            end

            % Generate an offspring
            Offspring = Global.Variation(Population(P(1:2)),1);

            % Update the ideal point
            Z = min(Z,Offspring.obj);

            % Calculate the constraint violation of offspring and P
            CVO = sum(max(0,Offspring.con));
            CVP = sum(max(0,Population(P).cons),2);
            
            % Update the solutions in P by PBI approach
            normW   = sqrt(sum(W(P,:).^2,2));
            normP   = sqrt(sum((Population(P).objs-repmat(Z,length(P),1)).^2,2));
            normO   = sqrt(sum((Offspring.obj-Z).^2,2));
            CosineP = sum((Population(P).objs-repmat(Z,length(P),1)).*W(P,:),2)./normW./normP;
            CosineO = sum(repmat(Offspring.obj-Z,length(P),1).*W(P,:),2)./normW./normO;
            g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
            g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
            Population(P(find(g_old>=g_new & CVP==CVO | CVP>CVO,nr))) = Offspring;
        end
    end
end