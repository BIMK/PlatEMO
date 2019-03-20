function hpaEA(Global)
% <algorithm> <H>
% Hyperplane assisted evolutionary algorithm

%------------------------------- Reference --------------------------------
% H. Chen, Y. Tian, W. Pedrycz, G. Wu, R. Wang, and L. Wang, Hyperplane
% assisted evolutionary algorithm for many-objective optimization problems,
% IEEE Transactions on Cybernetics, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen

    % Randomly initiallize a population
    Population = Global.Initialization();    
    % Generate a set of reference vectors
    [V, ~] = UniformPoint(Global.N, Global.M);
    % Record the nadir point
    maxObj = Inf(1, Global.M);
    % Initialize an array to record the indexes of prominent solutions
    PSI = [];
    
    %% Optimization
    while Global.NotTermination(Population)
        
        % Mating selection
        domNum = size(PSI, 2);
        MatingPool = randi(numel(Population), 1, Global.N-domNum);
        MatingPool = [MatingPool, PSI];
        MatingPool(randperm(size(MatingPool, 2))) = MatingPool(1: end);
        
        Offspring  = GA(Population(MatingPool));
        
        
        CombinePop = [Population, Offspring];
        % Remove some solutions that cannot dominate the worst point
        fitness = max(CombinePop.objs - repmat(maxObj, numel(CombinePop), 1), [], 2);
        deleteI = fitness > 0;
        CombinePop(deleteI) = [];
        if numel(CombinePop) < 2
            CombinePop = [Population, Offspring];
        end
        % Update the worst point
        tempMax = max(CombinePop.objs);
        replaceI = maxObj > tempMax;
        maxObj(replaceI) = tempMax(replaceI);        
        
        % Environmental Selection
        [Population, PSI] = hpaEnvironmentalSelection(CombinePop, Global, V);
          
    end
end