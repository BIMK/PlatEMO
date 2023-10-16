function Population = archive(Population,N)
% Select feasible and non-dominated solutions by using NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select feasible solutions
    fIndex           = all(Population.cons <= 0,2);
    Population       = Population(fIndex);

    if isempty(Population)
        return
    else
        %% Non-dominated sorting
        [FrontNo,~] = NDSort(Population.objs,1);
        Next = (FrontNo == 1);    
        Population = Population(Next);    
        if sum(Next) > N
            %% Calculate the crowding distance of each solution
            CrowdDis = CrowdingDistance(Population.objs);
            [~,Rank] = sort(CrowdDis,'descend');
            Population = Population(Rank(1:N));
        end

    end
end

% Calculate the crowding distance of each solution
function CrowdDis = CrowdingDistance(PopObj)
    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);

    for i = 1 : M
        [~,Rank] = sortrows(PopObj(:,i));
        CrowdDis(Rank(1))   = inf;
        CrowdDis(Rank(end)) = inf;

        for j = 2 : (N - 1)
            CrowdDis(Rank(j)) = CrowdDis(Rank(j))+(PopObj(Rank(j+1),i)-PopObj(Rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end