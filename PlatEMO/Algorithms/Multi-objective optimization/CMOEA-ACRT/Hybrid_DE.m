function [Offspring,M_recod] = Hybrid_DE(Population,Problem,MP_pool)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    [N,D]        = size(Population(1).decs);    
    [FrontNo1,~] = NDSort(Population.objs,Population.cons,inf);
    CrowdDis1    = CrowdingDistance(Population.objs,FrontNo1);
  
    [~,r1]    = sortrows([FrontNo1',-CrowdDis1']);
    R_sum(r1) = 1 : Problem.N;  % obtain ranking based on CDP
    [~,Rank]  = sort(R_sum);
    trial     = zeros(1*Problem.N,D);
    
    M_recod = zeros(Problem.N,1);
    for i = 1 : Problem.N
        rd      = rand;
        pro_sel = MP_pool(randi(5));
        
        if  rd <= pro_sel
            M_recod(i) = 1;
            cand_pool  = [0.1,0.2,1.0];
            CR         = cand_pool(randi(3));
            num_best   = 3;
            rd_index   = randperm(Problem.N,num_best);
            rd_index   = sort(rd_index);
            xr1        = Rank(rd_index(1));
            xr2        = Rank(rd_index(2));
            xr3        = Rank(rd_index(3));

            v      = Population(xr1).decs+rand*(Population(xr2).decs-Population(xr3).decs);%Mutation operator I
            Lower  = repmat(Problem.lower,N,1);
            Upper  = repmat(Problem.upper,N,1);
            v      = min(max(v,Lower),Upper);
            % crossover operation
            Site   = rand(N,D) < CR;
            j_rand = floor(rand * D) + 1;
            Site(1, j_rand) = 1;
            Site_  = 1-Site;
            trial(i, :)     = Site.*v+Site_.*Population(i).decs;
        else
            M_recod(i) = 2;
            cand_pool  = [0.1,0.8,1.0];
            F          = cand_pool(randi(3));

            indexset    = 1:Problem.N;
            indexset(i) = [];
            r1  = floor(rand*(Problem.N-1))+1;
            xr1 = indexset(r1);
            indexset(r1) = [];
            r2  = floor(rand*(Problem.N-2))+1;
            xr2 = indexset(r2);
            indexset(r2) = [];
            r3  = floor(rand*(Problem.N-3))+1;
            xr3 = indexset(r3);

            v      = Population(i).decs+F*(Population(xr1).decs-Population(i).decs)+F*(Population(xr2).decs-Population(xr3).decs);%Mutation operator II
            Lower  = repmat(Problem.lower,N,1);
            Upper  = repmat(Problem.upper,N,1);
            v      = min(max(v,Lower),Upper);
            trial(i, :) = v;
        end
    end 
    Offspring = trial;
    Offspring = Problem.Evaluation(Offspring);
end