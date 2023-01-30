classdef MaOEAIT < ALGORITHM
% <multi/many> <real/integer> <constrained/none>
% Many-objective evolutionary algorithms based on an independent two-stage
% approach
% Evaluation1 --- 20000 --- Number of evaluations for NDWA
% Evaluation1 ---  6000 --- Number of evaluations for reference lines mapping
% epsilon     --- 0.999 --- Threshold in PCA

%------------------------------- Reference --------------------------------
% Y. Sun, B. Xue, M. Zhang, G. G. Yen, A new two-stage evolutionary
% algorithm for many-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2019, 23(5): 748-761.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [Evaluation1,Evaluation2,epsilon] = Algorithm.ParameterSet(20000,6000,0.999);

            %% Generate random population
            Population = Problem.Initialization();
            Archive    = [];

            %% Generate the weight vectors
            [W,W_num] = UniformPoint(Problem.N,Problem.M);
            UW        = [W;W(end:-1:1,:)];
            num_W     = size(UW,1);

            %% The objective value of each solution on single-objective in different stage
            SF = zeros(Problem.N+floor(Problem.N/2)*2,1);

            %% The objective value of each solution on single-objective in phase 1
            SF(1:Problem.N) = sum(Population.objs.*UW(1,:),2);

            %% Non-dominated dynamic weight aggregation (NDWA)
            while Algorithm.NotTerminated(Population) && Problem.FE < Evaluation1
                W_index = mod(ceil(Problem.FE/Problem.N)+1,num_W);
                if W_index == 0
                    W_index = num_W;
                end
                MatingPool          = TournamentSelection(2,Problem.N,SF(1:Problem.N));
                Offspring           = OperatorGA(Problem,Population(MatingPool),{0.9,20,1,20});
                SF(Problem.N+1:end) = sum(Offspring.objs.*UW(W_index,:),2);
                [Population,SF]     = NDWA_EnvironmentalSelection([Population,Offspring],SF,Problem.N);
                Archive             = [Archive,Population];
                Archive             = Archive(NDSort(Archive.objs,1)==1);
            end

            %% Pareto-optimal subspace learning
            learned_data = FindSubspace(Archive.decs,epsilon);
            mean_val     = mean(learned_data,1);
            u_limit      = ones(size(Problem.upper));
            l_limit      = zeros(size(Problem.lower));
            for i = 1:numel(mean_val)
                if abs(learned_data(1,i)-mean_val(i)) < 0.1
                    l_limit(i) = round(mean_val(i),1);
                    u_limit(i) = round(mean_val(i),1);
                else
                    l_limit(i) = Problem.lower(i);
                    u_limit(i) = Problem.upper(i);
                end
            end

            %% Reference lines mapping
            % Optimize Problem.M single-objective optimization problems
            rf_list       = diag(ones(1,Problem.M));
            SO_maxeval    = floor(Evaluation2/Problem.M);
            Extreme_point = repmat(Population(1),1,Problem.M);
            for i = 1 : Problem.M
                Curr_evl        = Problem.FE;
                PopDec          = rand(Problem.N, Problem.D).*repmat(u_limit-l_limit,Problem.N,1) + repmat(l_limit,Problem.N,1);
                Population      = Problem.Evaluation(PopDec);
                SF(1:Problem.N) = cos_v_func(Population.objs,rf_list(i,:));
                while Algorithm.NotTerminated(Population) && Problem.FE < (Curr_evl+SO_maxeval)
                    MatingPool          = TournamentSelection(2,Problem.N,SF(1:Problem.N));
                    Offspring           = Operator(Problem,Population(MatingPool).decs,l_limit,u_limit);            
                    SF(Problem.N+1:end) = cos_v_func(Offspring.objs,rf_list(i,:));
                    [~,Rank]            = sort(SF,1);
                    Population          = [Population,Offspring];
                    Population          = Population(Rank(1:Problem.N));
                    SF(1:Problem.N)     = SF(Rank(1:Problem.N));
                end
                Extreme_point(i) = Population(1);
            end
            Ideal_point = min(Extreme_point.objs,[],1);
            Nadir_point = max(Extreme_point.objs,[],1);
            RefPoint    = W.*repmat(Nadir_point-Ideal_point,size(W,1),1) + repmat(Ideal_point,size(W,1),1);

            %% Diversity maintaining
            % Optimize W_num single-objective optimization problems
            SO_maxeval = floor((Problem.maxFE-Evaluation1-Evaluation2)/W_num);
            Result     = repmat(Population(1),1,size(RefPoint,1));
            for i = 1 : W_num
                Curr_evl        = Problem.FE;
                PopDec          = rand(Problem.N,Problem.D).*repmat(u_limit-l_limit,Problem.N, 1) + repmat(l_limit,Problem.N,1);
                Population      = Problem.Evaluation(PopDec);
                SF(1:Problem.N) = cos_v_func(Population.objs,RefPoint(i,:));
                while Algorithm.NotTerminated(Population) && Problem.FE < (Curr_evl+SO_maxeval)
                    MatingPool          = TournamentSelection(2,Problem.N,SF(1:Problem.N));
                    Offspring           = Operator(Problem,Population(MatingPool).decs,l_limit,u_limit);
                    SF(Problem.N+1:end) = cos_v_func(Offspring.objs,RefPoint(i,:));
                    [~,Rank]            = sort(SF,1);
                    Population          = [Population,Offspring];
                    Population          = Population(Rank(1:Problem.N));
                    SF(1:Problem.N)     = SF(Rank(1:Problem.N));
                    if i == W_num && Problem.FE == (Curr_evl+SO_maxeval)  
                        Result(W_num) = Population(1);
                        Population    = Result;
                    end 
                end
                Result(i) = Population(1);
            end
        end
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

function Offspring = Operator(Problem,ParentDec,lower,upper)
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

    Offspring = Problem.Evaluation(Offspring);
end