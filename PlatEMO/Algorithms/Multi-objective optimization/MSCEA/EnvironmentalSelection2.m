function [Population, Fitness] = EnvironmentalSelection2(Population,N,epsn)
% The environmental selection of MSCEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    [~, nCon] = size(Population.cons);
    PopCon    = max(0,Population.cons);
    if sum(sum(PopCon<=epsn, 2)==nCon) > N       
        tmp        = sum(PopCon<=epsn, 2)==nCon;
        Population = Population(1:end, tmp);
        CV         = sum(max(0,Population.cons),2);        
        Fitness    = CalFitness([Population.objs,CV]);        
    else       
        CV         = sum(max(0,Population.cons),2);
        Fitness    = CalFitness([CalSDE(Population.objs)',CV]);        
    end
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).objs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    Population = Population(Next);
    Fitness    = Fitness(Next);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end