function [Population,Fitness] = EnvironmentalSelectionSup(Population,N,net)
% The environmental selection of sup population using net

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    %% Calculate max min distance between reference points in V
    V = net.w;
    C = net.C;
    Distance = pdist2(V,V);
    Distance = Distance.*~C;
    d = zeros(1,length(V));
    for i = 1:length(V)
        Dis = Distance(i,:);
        Dis(Dis==0) = [];
        d(i) = min(Dis);
    end
    theta = max(d);
    
    %% Calculate the fitness of each solution
    Fitness = CalFitnessSup(Population.decs,V);
    
    %% Environmental selection
    Next = Fitness < theta;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).decs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
    % Sort the population
    [Fitness,rank] = sort(Fitness);
    Population = Population(rank);
end

function Del = Truncation(PopDec,K)
% Select part of the solutions by truncation

    %% Truncation
    D = pdist2(PopDec,PopDec);
    D(logical(eye(length(D)))) = inf;
    
    Del = false(1,size(PopDec,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(D(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end