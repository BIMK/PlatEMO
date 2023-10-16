function [Offspring] = DEgenerator2(Problem,Population)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    cv = sum(max(0,Population.cons),2);       

    FrontNo = NDSort(Population.objs,Population.cons,1);   
    index1  = find(FrontNo==1);
    r       = floor(rand*length(index1))+1;
    best    = index1(r);

    [N,D] = size(Population(1).decs);       
    trial = zeros(1*Problem.N,D);
       
    for i = 1 : Problem.N   
        % DE/rand-to-bestf/1/bin--Convergence
        if rand > 0.5
            l = rand;
            if l <= 1/3 
            	F = .6;
            elseif l <= 2/3
                F = 0.8;
            else
                F = 1.0;
            end   
            l = rand;
            if l <= 1/3
                CR = .1;
            elseif l <= 2/3
                CR = 0.2;
            else
                CR = 1.0;
            end
            indexset    = 1 : Problem.N;
            indexset(i) = [];
            r1  = floor(rand*(Problem.N-1))+1;
            xr1 = indexset(r1);
            indexset(r1) = [];
            r2  = floor(rand*(Problem.N-2))+1;
            xr2 = indexset(r2)  ;
            r3  = floor(rand*(Problem.N-3))+1;
            xr3 = indexset(r3);
            Best_index = Population(best).decs;
            v      = Population(xr1).decs+rand*(Best_index-Population(xr1).decs)+F*(Population(xr2).decs-Population(xr3).decs);  
            Lower  = repmat(Problem.lower,N,1);
            Upper  = repmat(Problem.upper,N,1);
            v      = min(max(v,Lower),Upper);
            Site   = rand(N,D) < CR;
            j_rand = floor(rand * D) + 1;
            Site(1, j_rand) = 1;
            Site_  = 1-Site;
            trial(i, :) = Site.*v+Site_.*Population(i).decs;  
        % DE/current-to-rand/1--Diversity
        else
            l = rand;
            if l <= 1/3
                F = .6;
            elseif l <= 2/3
                F = 0.8;
            else
                F = 1.0;
            end
            indexset    = 1:Problem.N;
            indexset(i) = [];
            r1  = floor(rand*(Problem.N-1))+1;
            xr1 = indexset(r1);
            indexset(r1) = [];
            r2  = floor(rand*(Problem.N-2))+1;
            xr2 = indexset(r2);
            indexset(r2) = [];
            r3    = floor(rand*(Problem.N-3))+1;
            xr3   = indexset(r3);
            v     = Population(i).decs+rand*(Population(xr1).decs-Population(i).decs)+F*(Population(xr2).decs-Population(xr3).decs); 
            Lower = repmat(Problem.lower,N,1);
            Upper = repmat(Problem.upper,N,1); 
            trial(i, :) = min(max(v,Lower),Upper);   
        end
    end
    Offspring = trial;
    Offspring = Problem.Evaluation(Offspring);
end