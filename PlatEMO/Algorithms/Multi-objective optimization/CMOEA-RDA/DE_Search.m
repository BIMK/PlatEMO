function Offspring= DE_Search(Population,Problem,pop_num)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    pop_num_total = length(Population);
    FrontNo       = NDSort(Population.objs,0,1);
    
    CrowdDis  = CrowdingDistance(Population.objs,FrontNo);
    [~,r1]    = sortrows([FrontNo',-CrowdDis']);
    R_sum(r1) = 1 : pop_num_total;  % obtain ranking based on CDP
    [~,Rank]  = sort(R_sum);
    
    [N,D] = size(Population(1).decs);
    trial = zeros(1*pop_num,D);
    
    for i = 1 : pop_num
        % Randomly select one individual from ranking level 1 as the best
        index_No1  = find(FrontNo==1);
        r          = floor(rand*length(index_No1))+1;
        best_index = index_No1(r);
    
        % Randomly select three individuals who are not equal to each other
        indexset    = 1 : pop_num_total;
        indexset(i) = [];
        r1  = floor(rand*(pop_num_total-1))+1;
        xr1 = indexset(r1);
        indexset(r1) = [];
        r2  = floor(rand*(pop_num_total-2))+1;
        xr2 = indexset(r2);
        indexset(r2) = [];
        r3  = floor(rand*(pop_num_total-3))+1;
        xr3 = indexset(r3);
    
        F_pool = [0.1,0.8,1.0];
        F      = F_pool(randi(3));
        CR     = 0.01;
    
        if rand <= 2/3
            if rand < 0.5
                v = Population(xr1).decs+F*(Population(best_index).decs-Population(xr1).decs)+F*(Population(xr2).decs-Population(xr3).decs);%Mutation operation 1
            else
                v = Population(xr1).decs+F*(Population(i).decs-Population(xr1).decs)+F*(Population(xr2).decs-Population(xr3).decs);%Mutation operation 2
            end
    
            % Boundary Repair
            Lower = repmat(Problem.lower,N,1);
            Upper = repmat(Problem.upper,N,1);
            v     = min(max(v,Lower),Upper);
    
            % Crossover operation
            Site   = rand(N,D) < CR;
            j_rand = floor(rand * D) + 1;
            Site(1, j_rand) = 1;
            Site_           = 1-Site;
            trial(i, :)     = Site.*v+Site_.*Population(i).decs;
        else
            num_best = 3;
            rd_index = randperm(pop_num_total,num_best);
            rd_index = sort(rd_index);
            xr1      = Rank(rd_index(1));
            xr2      = Rank(rd_index(2));
            xr3      = Rank(rd_index(3));
    
            v     = Population(xr1).decs+rand*(Population(xr2).decs-Population(xr3).decs);
            Lower = repmat(Problem.lower,N,1);
            Upper = repmat(Problem.upper,N,1);
            trial(i, :) = min(max(v,Lower),Upper);
        end
    end
    Offspring = trial;
    Offspring = Problem.Evaluation(Offspring);
end