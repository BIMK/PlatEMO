function [ Offspring ] = CMMS(Population, Problem, p,Population2, flag_index,Direc_index)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    Nw          = min(Problem.D,Problem.N);
    permutation = randperm(Problem.N);
    Population  = Population(permutation);
    
    N = ceil(Problem.N*p);
    D = Problem.D;
    
    for i = 1 : N
        r1(i) = randi([1,Problem.N],1);
        while r1(i) == i
            r1(i) = randi([1,Problem.N],1);
        end
    
        r2(i) = randi([1,Problem.N],1);
        while r2(i)==i || r2(i)==r1(i)
            r2(i) = randi([1,Problem.N],1);
        end
    
        r3(i) = randi([1,Problem.N],1);
        while r3(i)==r2(i) || r3(i)==r1(i) || r3(i)==i
            r3(i) = randi([1,Problem.N],1);
        end
    end
    Offspring = [];
    for j = Direc_index
        start_X = Population(1:N).decs;
        if flag_index == 2
            end_X = Population2(r1).decs;
            if j == 2
                end_X = Problem.lower+Problem.upper - end_X;
            end
        elseif flag_index == 1
            end_X = Population(r2).decs;
            if j == 2
                end_X = Problem.lower+Problem.upper - end_X;
            end
        end
        Intervalmax = sqrt(sum((Problem.upper-Problem.lower).^2,2));
        Intervalmin = 0;
        RandSample  = Intervalmin + rand(N,Nw)*(Intervalmax-Intervalmin);
    
        for i = 1 : Nw
            result_pop = end_X + repmat(RandSample(:,i),1,D).* (start_X - end_X);
            result_pop = boundConstraint (result_pop, Population(1:N).decs, [Problem.lower;Problem.upper]);
            Offspring  = [Offspring;result_pop];
        end
    end
    %% Polynomial mutation
    [proM,disM] = deal(1,20);
    Lower       = repmat(Problem.lower,N,1);
    Upper       = repmat(Problem.upper,N,1);
    Site        = rand(N,D) < proM/D;
    mu          = rand(N,D);
    temp        = Site & mu<=0.5;
    
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    
    Offspring = Problem.Evaluation(Offspring);
end