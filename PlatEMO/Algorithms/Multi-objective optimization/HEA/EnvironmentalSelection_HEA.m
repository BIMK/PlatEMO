function [NewPopulation, Population_hd] = EnvironmentalSelection_HEA(Population, hd, zmin, zmax, N, W)
% The environmental selection of HEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu

    %% Normalization the solutions
    PopObj = Population.objs;
    [P, ~] = size(PopObj);
    zmin       = min(zmin, min(PopObj,[],1));
    nf         = zmax - zmin;
    % Avoid dividing by zero
    nf(nf < 1e-6)   = 1e-6;                     
    PopObj       = (PopObj - zmin) ./ nf;
    
    %% Environmental selection
    Fitness = pdist2(PopObj,W,'cosine');
    [~, index] = min(Fitness, [], 1);
    Chosen = zeros(P, 1);
    for i = 1: N
        Chosen(index(i)) = 1;
    end
    NewPopulation = Population(Chosen == 1);
    Population_hd = hd(Chosen == 1);
    
    %% Solutions supplement operation
    K = N - sum(Chosen);
    for i = 1: K
        [~, Max_hd_index] = max(hd);
        NewPopulation = [NewPopulation, Population(Max_hd_index)];
        Population_hd = [Population_hd, hd(Max_hd_index)];
        hd(Max_hd_index) = -inf;
    end
end