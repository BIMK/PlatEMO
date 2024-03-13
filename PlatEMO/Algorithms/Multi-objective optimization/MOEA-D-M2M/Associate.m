function Population = Associate(Population,W,S)
% Allocation of solutions to subproblems

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    K = size(W,1);
    
    %% Allocation of solutions to subproblems
    % Transformation
    [~,transformation] = max(1-pdist2(Population.objs,W,'cosine'),[],2);
    partition = zeros(S,K);
    % Allocation
    for i = 1 : K
        current = find(transformation==i);
        if length(current) < S
            % Randomly select solutions and join to the current subproblem
            current = [current;randi(length(Population),S-length(current),1)];
        elseif length(current) > S
            % Delete solutions from the current subproblem by non-dominated
            % sorting and crowding distance
            [FrontNo,MaxFNo] = NDSort(Population(current).objs,S);
            Last = find(FrontNo==MaxFNo);
            CrowdDis = CrowdingDistance(Population(current(Last)).objs);
            [~,rank] = sort(CrowdDis);
            FrontNo(Last(rank(1:sum(FrontNo<=MaxFNo)-S))) = inf;
            current = current(FrontNo<=MaxFNo);
        end
        partition(:,i) = current;
    end
    Population = Population(partition(:));
end