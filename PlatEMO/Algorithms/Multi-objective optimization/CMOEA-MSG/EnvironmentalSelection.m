function [Population,Fitness] = EnvironmentalSelection(Population,N,index,count,c)
% Environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    a        = 0;
    Fitness1 = [];
    num      = 1;
    Pop      = [];
    PopObj   = Population.objs;
	Fitness  = CalFitness(Population.objs,Population.cons,index,count,c);
    M        = size(Population(1).objs,2);
    L        = size(Population,2);
    rank     = inf(1,L);
    
    z          = min(PopObj,[],1);
    [W,~]      = UniformPoint(100,M);
    [~,Region] = min(pdist2(PopObj-z,W,'cosine'),[],2);  
    for i = 1 : size(W,1)
        index = find(Region==i);
         if ~isempty(index)
            [~,index1] = sort(Fitness(index));
            index = index(index1);
            for j = 1 : size(index)
                rank(index(j)) = j;
            end
        end
    end
    while true
        [~,index3] = find(rank==num);
        a = a + length(index3);
        if a <= N
            Pop = [Pop,Population(index3)];
            Fitness1 = [Fitness1,Fitness(index3)];
            num = num + 1;
        else
            break;
        end
    end
    if length(Pop) == N
        Population = Pop;
        Fitness    = Fitness1;
    else
        [~,index3] = find(rank==num);
        Pop1       = Population(index3);
        Fitness2   = Fitness(index3);
        Del        = Truncation(Pop1.objs,a-N);
        Population = [Pop,Pop1(~Del)];
        Fitness    = [Fitness1,Fitness2(~Del)];
    end
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

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