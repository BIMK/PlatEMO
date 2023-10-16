function [Population,FrontNo,DWeight] = EnvironmentalSelection(Population,N)
% The environmental selection of DWU

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Gladston Moreira

    %% Non-dominated sorting
    [FrontNo] = NDSort(Population.objs,Population.cons,N);
    
    %% Calculate the dominance information each solution
    DWeight = InfoDominance(Population.objs);

    %% Environment Select
    Next = ReplacementUniformity(Population,N,FrontNo,DWeight);
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    DWeight    = DWeight(Next);
end

function InfoD = InfoDominance(PopObj)
% Calculate the information dominance each solution

    N = size(PopObj,1);

    %% Dominance count each solution
    D = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
            if k == 1
                D(i,j) = true;
            elseif k == -1
                D(j,i) = true;
            end
        end
    end
    CountDominance = sum(D,2);
    
    %% Calculate information dominance each solution
    InfoD = D'*CountDominance;
end