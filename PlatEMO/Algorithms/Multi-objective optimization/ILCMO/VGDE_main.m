function Offspring = VGDE_main(Problem, Population, Fitness,numberOfGroups, typeOfGroups,P)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N         = length(Population);
    Offspring = [];
    
    %% Mating selection
    MatingPool = TournamentSelection(2,N,Fitness);
    Population = Population(MatingPool);
    Fitness    = Fitness(MatingPool);
    
    %% Clustering
    [~, indBest] = sort(Fitness, 'descend');
    S1           = Population(indBest(1:N/2));%Winner
    location1    = rand(1,N/2)< P;
    S2           = Population(indBest(1+N/2:end));%Loser
    location2    = rand(1,N/2)<= P;
    
    %% Diversity evolution
    r0         = randperm(N/2);
    [r1,r2,r3] = gnR1R2R3(N/2, r0);
    Parent11   = S1(r0);
    Parent12   = S1(r1);
    Parent13   = S1(r2);
    
    if ~isempty(find(location1==1, 1))
        Off       = Group_DE_rand_1(Problem,Parent11(location1).decs,Parent12(location1).decs,Parent13(location1).decs,numberOfGroups, typeOfGroups);
        Offspring = [Offspring,Off];
    end
    if ~isempty(find(location1==0, 1))
        Off       = DE_rand_1(Problem,Parent11(~location1).decs,Parent12(~location1).decs,Parent13(~location1).decs);
        Offspring = [Offspring,Off];
    end

    %% Convergence evolution 
    Parent21 = S2(randperm(N/2));
    Parent22 = S1(randperm(N/2));
    if ~isempty(find(location2==1, 1))
        Off       = Group_DE_best_1(Problem,Parent21(location2).decs,Parent22(location2).decs,numberOfGroups, typeOfGroups);
        Offspring = [Offspring,Off];
    end
    if ~isempty(find(location2==0, 1))
        Off       = DE_best_1(Problem,Parent21(~location2).decs,Parent22(~location2).decs);
        Offspring = [Offspring,Off];
    end
end

function Offspring = Group_DE_rand_1(Problem,Parent1,Parent2,Parent3,numberOfGroups, typeOfGroups)
    [N,D] = size(Parent1);
    Fm    = [0.6,0.8,1.0];
    index = randi([1,length(Fm)],N,1);
    F     = Fm(index);
    F     = F';
    F     = F(:, ones(1,D));

    %% Differental evolution
    [outIndexList,~] = CreateGroups(numberOfGroups,Parent1,Problem.D,typeOfGroups);
    chosengroups     = randi(numberOfGroups,size(outIndexList,1),1);
    Site             = outIndexList == chosengroups;
    Offspring        = Parent1;
    Offspring(Site)  = Offspring(Site) + F(Site).*(Parent2(Site)-Parent3(Site));
    
    [proM,disM] = deal(1,20);
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    mu    = rand(N,1);
    mu    = repmat(mu,1,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = Problem.Evaluation(Offspring);
end

function Offspring = DE_rand_1(Problem,Parent1,Parent2,Parent3)
    [N,D] = size(Parent1);
    Fm    = [0.6,0.8,1.0];
    CRm   = [0.1,0.2,1.0];
    index = randi([1,length(Fm)],N,1);
    F     = Fm(index);
    F     = F';
    F     = F(:, ones(1,D));
    index = randi([1,length(CRm)],N,1);
    CR    = CRm(index);
    CR    = CR';

    %% Differental evolution
    Site            = rand(N,D) < repmat(CR,1,D);
    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + F(Site).*(Parent2(Site)-Parent3(Site));

    %% Polynomial mutation
    [proM,disM] = deal(1,20);
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = Problem.Evaluation(Offspring);
end

function Offspring = Group_DE_best_1(Problem,Parent1,Parent2,numberOfGroups,typeOfGroups)
    [N,D] = size(Parent1);
    Fm    = [0.6,0.8,1.0];
    index = randi([1,length(Fm)],N,1);
    F     = Fm(index);
    F     = F';
    F     = F(:, ones(1,D));

    %% Differental evolution
    [outIndexList,~] = CreateGroups(numberOfGroups,Parent1,D,typeOfGroups);
    chosengroups     = randi(numberOfGroups,size(outIndexList,1),1);
    Site             = outIndexList == chosengroups;
    Offspring        = Parent1;
    Offspring(Site)  = Offspring(Site) + F(Site).*(Parent2(Site)-Offspring(Site));

    %% Polynomial mutation
    [proM,disM] = deal(1,20);
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    mu    = rand(N,1);
    mu = repmat(mu,1,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = Problem.Evaluation(Offspring);
end

function Offspring = DE_best_1(Problem,Parent1,Parent2)
    [N,D] = size(Parent1);
    Fm    = [0.6,0.8,1.0];
    CRm   = [0.1,0.2,1.0];
    index = randi([1,length(Fm)],N,1);
    F     = Fm(index);
    F     = F';
    F     = F(:, ones(1,D));
    index = randi([1,length(CRm)],N,1);
    CR    = CRm(index);
    CR    = CR';

    %% Differental evolution
    Site            = rand(N,D) < repmat(CR,1,D);
    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + F(Site).*(Parent2(Site)-Offspring(Site));

    %% Polynomial mutation
    [proM,disM] = deal(1,20);
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = Problem.Evaluation(Offspring);
end