function MaOEAIT(Global)
% <algorithm> <M>
% Many-objective evolutionary algorithms based on an independent two-stage
% approach
% Evaluation1 --- 20000 --- Number of evaluations for NDWA
% Evaluation1 ---  6000 --- Number of evaluations for reference lines mapping
% epsilon     --- 0.999 --- Threshold in PCA

%------------------------------- Reference --------------------------------
% Y. Sun, B. Xue, M. Zhang, G. G. Yen, A new two-stage evolutionary
% algorithm for many-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2018.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [Evaluation1,Evaluation2,epsilon] = Global.ParameterSet(20000,6000,0.999);

    %% Generate random population
    Population = Global.Initialization();
    Archive    = [];
    
    %% Generate the weight vectors
    [W,W_num] = UniformPoint(Global.N,Global.M);
    UW        = [W;W(end:-1:1,:)];
    num_W     = size(UW,1);

    %% The objective value of each solution on single-objective in different stage
    SF = zeros(Global.N+floor(Global.N/2)*2,1);

    %% The objective value of each solution on single-objective in phase 1
    SF(1:Global.N) = sum(Population.objs.*UW(1,:),2);

    %% Non-dominated dynamic weight aggregation (NDWA)
    while Global.NotTermination(Population) && Global.evaluated < Evaluation1
        W_index = mod(Global.gen+1,num_W);
        if W_index == 0
            W_index = num_W;
        end
        MatingPool         = TournamentSelection(2,Global.N,SF(1:Global.N));
        Offspring          = GA(Population(MatingPool),{0.9,20,1,20});
        SF(Global.N+1:end) = sum(Offspring.objs.*UW(W_index,:),2);
        [Population,SF]    = NDWA_EnvironmentalSelection([Population,Offspring],SF,Global.N);
        Archive            = [Archive,Population];
        Archive            = Archive(NDSort(Archive.objs,1)==1);
    end

    %% Pareto-optimal subspace learning
    learned_data = FindSubspace(Archive.decs,epsilon);
    mean_val     = mean(learned_data,1);
    u_limit      = ones(size(Global.upper));
    l_limit      = zeros(size(Global.lower));
    for i = 1:numel(mean_val)
        if abs(learned_data(1,i)-mean_val(i)) < 0.1
            l_limit(i) = round(mean_val(i),1);
            u_limit(i) = round(mean_val(i),1);
        else
            l_limit(i) = Global.lower(i);
            u_limit(i) = Global.upper(i);
        end
    end

    %% Reference lines mapping
    % Optimize Global.M single-objective optimization problems
    rf_list       = diag(ones(1,Global.M));
    SO_maxeval    = floor(Evaluation2/Global.M);
    Extreme_point = repmat(Population(1),1,Global.M);
    for i = 1 : Global.M
        Curr_evl       = Global.evaluated;
        PopDec         = rand(Global.N, Global.D).*repmat(u_limit-l_limit,Global.N,1) + repmat(l_limit,Global.N,1);
        Population     = INDIVIDUAL(PopDec);
        SF(1:Global.N) = cos_v_func(Population.objs,rf_list(i,:));
        while Global.NotTermination(Population) && Global.evaluated < (Curr_evl+SO_maxeval)
            MatingPool         = TournamentSelection(2,Global.N,SF(1:Global.N));
            Offspring          = EA_real(Population(MatingPool).decs,l_limit,u_limit);            
            SF(Global.N+1:end) = cos_v_func(Offspring.objs,rf_list(i,:));
            [~,Rank]           = sort(SF,1);
            Population         = [Population,Offspring];
            Population         = Population(Rank(1:Global.N));
            SF(1:Global.N)     = SF(Rank(1:Global.N));
        end
        Extreme_point(i) = Population(1);
    end
    Ideal_point = min(Extreme_point.objs,[],1);
    Nadir_point = max(Extreme_point.objs,[],1);
    RefPoint    = W.*repmat(Nadir_point-Ideal_point,size(W,1),1) + repmat(Ideal_point,size(W,1),1);
    
    %% Diversity maintaining
    % Optimize W_num single-objective optimization problems
    SO_maxeval = floor((Global.evaluation-Evaluation1-Evaluation2)/W_num);
    Result     = repmat(Population(1),1,size(RefPoint,1));
    for i = 1 : W_num
        Curr_evl       = Global.evaluated;
        PopDec         = rand(Global.N,Global.D).*repmat(u_limit-l_limit,Global.N, 1) + repmat(l_limit,Global.N,1);
        Population     = INDIVIDUAL(PopDec);
        SF(1:Global.N) = cos_v_func(Population.objs,RefPoint(i,:));
        while Global.NotTermination(Population) && Global.evaluated < (Curr_evl+SO_maxeval)
            MatingPool         = TournamentSelection(2,Global.N,SF(1:Global.N));
            Offspring          = EA_real(Population(MatingPool).decs,l_limit,u_limit);
            SF(Global.N+1:end) = cos_v_func(Offspring.objs,RefPoint(i,:));
            [~,Rank]           = sort(SF,1);
            Population         = [Population,Offspring];
            Population         = Population(Rank(1:Global.N));
            SF(1:Global.N)     = SF(Rank(1:Global.N));
            if i == W_num && Global.evaluated == (Curr_evl+SO_maxeval)  
                Result(W_num) = Population(1);
                Population    = Result;
            end 
        end
        Result(i) = Population(1);
    end
end

function [Population,SF] = NDWA_EnvironmentalSelection(Population,SF,N)
% Environmental selection

    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    Last = find(FrontNo==MaxFNo);
    Next(Last(1:N-sum(Next))) = true;
    Population = Population(Next);
    SF(1:N)    = SF(Next);
end

function SF = cos_v_func(PopObj,v)
% Calculate the objective value of each solution on each single-objective
% optimization problem

    r1 = sum(PopObj.*v,2);
    r2 = sqrt(sum(PopObj.^2,2));
    SF = -r1./r2;
end

function Offspring = EA_real(ParentDec,lower,upper)
% Simulated binary crossover and polynomial mutation

    [proC,disC,proM,disM] = deal(1,20,1,20);
    Parent1 = ParentDec(1:floor(end/2),:);
    Parent2 = ParentDec(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);

    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];

    Lower = repmat(lower,2*N,1);
    Upper = repmat(upper,2*N,1);
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

    Offspring = INDIVIDUAL(Offspring);
end