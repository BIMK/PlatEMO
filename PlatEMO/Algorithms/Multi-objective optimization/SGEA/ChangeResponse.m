function [Population,Archive,Centroid] = ChangeResponse(Problem,Population,Archive,Centroid)
% React to the environmental change

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select old solutions
    remain = ~Truncation(Population.objs,ceil(length(Population)/2));
    Population(remain) = Problem.Evaluation(Population(remain).decs);
    Ct = mean(Archive.decs,1);
    St = norm(Ct-Centroid);
    
    %% Generate new solutions
    Archive = Population(remain).best;
    CA = mean(Archive.decs,1);
    CR = mean(Population(remain).decs,1);
    X  = Population(~remain).decs;
    X  = X + repmat(St.*(CA-CR)./norm(CA-CR),size(X,1),1) + randn(size(X)).*St/2/sqrt(size(X,2));
    Population(~remain) = Problem.Evaluation(X);
    
    Archive  = Population.best;
    Centroid = Ct;
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