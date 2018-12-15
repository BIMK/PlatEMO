function Population = EnvironmentalSelection(Population,Weight,K)
% The environmental selection in MSOPS-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.objs;
    PopCon = Population.cons;

    %% Calculate the metric matrix
    [WMM,VADS] = CalMetric(PopObj,Weight);
    S = [WMM,VADS];
    
    %% Calculate the aggregate fitness of each solution
    % Scale each column by the minimum value found
    [minValue,minIndex] = min(S,[],1);
    S1 = S./repmat(minValue,size(S,1),1);
    % Scale the rows which have the minimum value by the second minimum
    % value
    for i = 1 : length(minIndex)
        S1(minIndex(i),i) = S(minIndex(i),i)./min(S([1:minIndex(i)-1,minIndex(i)+1:end],i));
    end
    % Calculate the aggregate fitness and sort the solutions
    r = min(S1,[],2);
    
    %% Identify the feasible solutions and the extreme solutions
    Feasible    = find(all(PopCon<=0,2));
    [~,extreme] = min(PopObj(Feasible,:),[],1);
    extreme     = Feasible(extreme);
    % Decrease the fitness value of extreme solutions so that they can be
    % selected preferentially
    r(extreme) = r(extreme) - (max(r)-min(r));
    % Set the fitness of each infeasible solution to its degree of
    % constraint violationm, then increase the value so that they will be
    % selected in the last
    Infeasible    = find(any(PopCon>0,2));
    r(Infeasible) = sum(PopCon(Infeasible,:).^2,2) + max(r);
    
    %% Select the best N solutions
    [~,rank] = sort(r);
    Next     = rank(1:K);
    % Population for next generation
    Population = Population(Next);
end