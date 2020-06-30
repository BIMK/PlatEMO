function MOEAPC(Global)
% <algorithm> <M>
% Multiobjective evolutionary algorithm based on polar coordinates
% delta --- 0.8 --- The probability of choosing parents locally
% T     ---  20 --- Neighborhood size 

%------------------------------- Reference --------------------------------
% R. Denysiuk, L. Costa, I. E. Santo, and J. C. Matos, MOEA/PC:
% Multiobjective evolutionary algorithm based on polar coordinates,
% Proceedings of the International Conference on Evolutionary
% Multi-Criterion Optimization, 2015, 141-155.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
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
            P = B(i,randperm(size(B,2),3));
        else
            P = randperm(nGrids,3);
        end  
        
        % Generate an offspring
        Offspring = DE(Population(P(1)),Population(P(2)),Population(P(3)));
        
        % Update the ideal point
        Z = min(Z,Offspring.obj);
        
        % Update the solutions in P by Environmental Selection
        Population = EnvironmentalSelection(Population,Offspring,Z,nDivs);
    end
end