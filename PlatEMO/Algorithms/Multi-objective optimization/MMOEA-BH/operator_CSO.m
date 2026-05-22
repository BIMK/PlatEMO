function Population = operator_CSO(Problem, pop_LBA, lower_bound, upper_bound, N, maxFEs)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    [V, N] = UniformPoint(N, Problem.M);
    
    Dec      = unifrnd(repmat(lower_bound, N, 1), repmat(upper_bound, N, 1));
    pop_init = Problem.Evaluation(Dec);
    
    if numel(pop_LBA) < N
        Population = [pop_LBA, pop_init];
        Population = Population(1:N);
    else
        Population = pop_LBA(1:N);
    end
    
    currentFE  = length(pop_init);
    Population = EnvironmentalSelection(Population, V, (currentFE/maxFEs)^2);
    
    %% Optimization
    while currentFE < maxFEs
        Fitness = calFitness(Population.objs);
        if length(Population) >= 2
            Rank = randperm(length(Population), floor(length(Population)/2)*2);
        else
            Rank = [1, 1];
        end
        
        Loser  = Rank(1:end/2);
        Winner = Rank(end/2+1:end);
        Change = Fitness(Loser) >= Fitness(Winner);
        Temp   = Winner(Change);
        
        Winner(Change) = Loser(Change);
        Loser(Change)  = Temp;
        
        % Generate offspring using the current subspace bounds
        Offspring = Operator_niche(Problem, Population(Loser), Population(Winner), lower_bound, upper_bound);
        
        currentFE  = currentFE + length(Offspring);
        Population = EnvironmentalSelection([Population, Offspring], V, (currentFE/maxFEs)^2);
    end
end

function Fitness = calFitness(PopObj)
% Calculate the fitness by shift-based density
    N      = size(PopObj,1);
    fmax   = max(PopObj,[],1);
    fmin   = min(PopObj,[],1);
    PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    Dis    = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    Fitness = min(Dis,[],2);
end

function Population = EnvironmentalSelection(Population,V,theta)
% The environmental selection of LMOCSO

    Population = Population(NDSort(Population.objs,1)==1);
    PopObj     = Population.objs;
    [N,M]      = size(PopObj);
    NV         = size(V,1);
    
    %% Translate the population
    PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
    
    %% Calculate the degree of violation of each solution
    CV     = sum(max(0,Population.cons),2);

    %% Calculate the smallest angle value between each vector and others
    cosine = 1 - pdist2(V,V,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2);

    %% Associate each solution to a reference vector
    Angle         = acos(1-pdist2(PopObj,V,'cosine'));
    [~,associate] = min(Angle,[],2);

    %% Select one solution for each reference vector
    Next = zeros(1,NV);
    for i = unique(associate)'
        current1 = find(associate==i & CV==0);
        current2 = find(associate==i & CV~=0);
        if ~isempty(current1)
            % Calculate the APD value of each solution
            APD      = (1+M*theta*Angle(current1,i)/gamma(i)).*sqrt(sum(PopObj(current1,:).^2,2));
            [~,best] = min(APD);
            Next(i)  = current1(best);
        elseif ~isempty(current2)
            % Select the one with the minimum CV value
            [~,best] = min(CV(current2));
            Next(i)  = current2(best);
        end
    end
    % Population for next generation
    Population = Population(Next(Next~=0));
end

function Offspring = Operator_niche(Problem, Loser, Winner, lower_bound, upper_bound)
% The competitive swarm optimizer of LMOCSO customized for local boundaries

    %% Parameter setting
    LoserDec  = Loser.decs;
    WinnerDec = Winner.decs;
    [N,D]     = size(LoserDec);
    LoserVel  = Loser.adds(zeros(N,D));
    WinnerVel = Winner.adds(zeros(N,D));
    
    %% Competitive swarm optimizer
    r1     = repmat(rand(N,1),1,D);
    r2     = repmat(rand(N,1),1,D);
    OffVel = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
    OffDec = LoserDec + OffVel + r1.*(OffVel-LoserVel);
    
    %% Add the winners
    OffDec = [OffDec;WinnerDec];
    OffVel = [OffVel;WinnerVel];
    
    %% Polynomial mutation (using the externally provided local block bounds)
    Lower = repmat(lower_bound, 2*N, 1);
    Upper = repmat(upper_bound, 2*N, 1);
    disM  = 20;
    Site  = rand(2*N,D) < 1/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    OffDec       = max(min(OffDec,Upper),Lower);
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
               
    temp         = Site & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = Problem.Evaluation(OffDec,OffVel);
end