function [Population,Fitness1] = KCCMO_sampling(Population,CA1,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    CA = CA1.best.objs;
    if isempty(CA)
        CA = CA1.objs;
    end

    Fitness = CalFitness_new(Population.obj,Population.con,CA);

    NCluster = N;
    [IDX,~] = kmeans(Population.obj,NCluster);
    Next = false(size(Population.obj,1),1);
    for i = 1:NCluster
        select = find(IDX == i);
        Fitness1 = Fitness(select);
        [~,index] = min(Fitness1);
        Next(select(index)) = true;
    end
    Population = givevalue(Population,Next);
end

function Fitness = CalFitness_new(PopObj,PopCon,CCA)

    N = size(PopObj,1);
    if isempty(PopCon)
        CV = zeros(N,1);
    else
        CV = sum(max(0,PopCon),2);
    end

    %% Detect the dominance relation between each two solutions
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            if CV(i) < CV(j)
                Dominate(i,j) = true;
            elseif CV(i) > CV(j)
                Dominate(j,i) = true;
            else
                k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
                if k == 1
                    Dominate(i,j) = true;
                elseif k == -1
                    Dominate(j,i) = true;
                end
            end

        end
    end

    S = sum(Dominate,2);

    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end

    Distance = pdist2(PopObj,CCA);
    Distance1 = min(Distance,[],2);
    D = 1./(Distance1+2);
    Fitness = R + D';
    end