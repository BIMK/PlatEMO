function Offspring = AES(Problem, Population, Fitness, Zmin)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
   
    N = length(Population);
    
    %% Mating selection
    MatingPool = TournamentSelection(2,N,Fitness);
    Population = Population(MatingPool);
    Fitness    = Fitness(MatingPool);
    
    %% Clustering
    [~,indBest] = sort(Fitness, 'ascend');
    Se          = Population(indBest(1:ceil(N/2)));     % Winner
    Sp          = Population(indBest(1+ceil(N/2):end)); % Loser

    % Perform bidirectional sampling on each solution in Se to generate new solutions.
    Off1 = DirectedSampling(Problem,Se);

    % For each solution in Sp, find the paired solution in Se that forms the minimum angle with it.
    [Mate1,Mate2] = Pair(Se,Sp,Zmin);
    Off2          = DE_best_1(Problem,Mate1.decs,Mate2.decs);

    Offspring = [Off1,Off2];
end

function [Mate1,Mate2] = Pair(Se,Sp,Zmin)

    %% Normalization
    SeObj   = Se.objs;
    [Num,M] = size(SeObj);
    SeObj   = (SeObj - repmat(Zmin,Num,1))./repmat(sqrt(sum(SeObj.^2,2)),1,M);
    SpObj   = Sp.objs;
    SpObj   = (SpObj - repmat(Zmin,Num,1))./repmat(sqrt(sum(SpObj.^2,2)),1,M);
    
    %% Association
    Cosine        = 1-pdist2(SpObj,SeObj,'cosine');
    [~,associate] = max(Cosine,[],2);
    Mate1         = Sp;
    Mate2         = Se(associate);
end

function Offspring = DE_best_1(Problem,Parent1,Parent2)
    [N,D] = size(Parent1);

    %% Differental evolution
    Intervalmin = sqrt(sum((Parent2-Parent1).^2,2));
    Intervalmax = sqrt(sum((Problem.upper-Problem.lower).^2,2));
    sigma       = Intervalmin + (Intervalmax-Intervalmin)*(1-Problem.FE/Problem.maxFE)^2;
    F           = rand(N,1).*sigma; F = F(:, ones(1,D));
    Offspring   = Parent1;
    Offspring   = Offspring + F.*(Parent2-Offspring)./Intervalmin;

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

function SampleSolution= DirectedSampling(Problem,Population)
    Nw     = length(Population);
    PopDec = Population.decs;
    Upper  = Problem.upper;
    Lower  = Problem.lower;

    Directnorm = [sqrt(sum((PopDec - repmat(Lower,Nw,1)).^2,2));sqrt(sum((PopDec - repmat(Upper,Nw,1)).^2,2))];
    Direction  = [PopDec - repmat(Lower,Nw,1);PopDec - repmat(Upper,Nw,1)]./repmat(Directnorm,1,Problem.D);

    %% Generate guiding solutions
    Intervalmax    = sqrt(sum((Upper-Lower).^2,2))*(1-Problem.FE/Problem.maxFE).^2;
    Intervalmin    = 0;
    RandSample     = Intervalmin + rand(1,Nw*4)*(Intervalmax-Intervalmin);
    SampleSolution = GenerateSampleSolution(Problem,PopDec,RandSample,Direction);

end

function SampleSolution = GenerateSampleSolution(Problem,X,RandSample,Direct)
% Generate some sample solutions along with the guiding directions

    Nw   = length(RandSample)/4;
    PopX = [X + repmat(RandSample(1:Nw)',1,Problem.D).* Direct(1:Nw,:);...
        X - repmat(RandSample(1+Nw:Nw*2)',1,Problem.D).* Direct(1:Nw,:);...
        X + repmat(RandSample(Nw*2+1:Nw*3)',1,Problem.D).* Direct(Nw+1:end,:);...
        X - repmat(RandSample(Nw*3+1:end)',1,Problem.D).* Direct(Nw+1:end,:)];
    PopX = max(min(repmat(Problem.upper,size(PopX,1),1),PopX),repmat(Problem.lower,size(PopX,1),1));
    SampleSolution = Problem.Evaluation(PopX);
end