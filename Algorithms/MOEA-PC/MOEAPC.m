function MOEAPC(Global)
% <algorithm> <H-N>
% MOEA/PC: Multiobjective Evolutionary Algorithm Based on Polar Coordinates
% delta --- 0.8 --- The probability of choosing parents locally
% T     ---  20 --- Neighborhood size 
% operator      --- DE

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

    %% Parameter setting
    [delta,T] = Global.ParameterSet(0.8,20);    

    %% Generate grid
    [G,nGrids,nDivs] = GridStructure(Global.N,Global.M);  
    Global.N = nGrids;

    %% Detect the neighbours of each grid
    B = pdist2(G,G);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    clear G

    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        % Select solution at random
        i = randi(nGrids);
        
        % Choose the parents
        if rand < delta
            P = B(i,randperm(size(B,2)));
        else
            P = randperm(nGrids);
        end
        parents = P(1:3);   
        
        % Generate an offspring
        Offspring = Global.Variation(Population(parents),1,@DE);
        
        % Update the ideal point
        Z = min(Z,Offspring.obj);
        
        % Update the solutions in P by Environmental Selection
        Population = EnvironmentalSelection(Population,Offspring,Z,nDivs);
    end
end