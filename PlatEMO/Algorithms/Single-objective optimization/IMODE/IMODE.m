classdef IMODE < ALGORITHM
% <single> <real/integer> <large/none> <constrained/none>
% Improved multi-operator differential evolution
% minN  ---   4 --- Minimum population size
% aRate --- 2.6 --- Ratio of archive size to population size

%------------------------------- Reference --------------------------------
% K. M. Sallam, S. M. Elsayed, R. K. Chakrabortty, and M. J. Ryan, Improved
% multi-operator differential evolution algorithm for solving unconstrained
% problems, Proceedings of the IEEE Congress on Evolutionary Computation,
% 2020.
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
            [minN,aRate] = Algorithm.ParameterSet(4,2.6);
            
            %% Generate random population
            Population = Problem.Initialization();
            Archive    = [];
            MCR = zeros(20*Problem.D,1) + 0.2;
            MF  = zeros(20*Problem.D,1) + 0.2;
            k   = 1;
            MOP = ones(1,3)/3;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Reduce the population size
                N          = ceil((minN-Problem.N)*Problem.FE/Problem.maxFE) + Problem.N;
                [~,rank]   = sort(FitnessSingle(Population));
                Population = Population(rank(1:N));
                Archive    = Archive(randperm(end,min(end,ceil(aRate*N))));
                % Generate parents, CR, F, and operator for each offspring
                Xp1 = Population(ceil(rand(1,N).*max(1,0.25*N))).decs;
                Xp2 = Population(ceil(rand(1,N).*max(2,0.5*N))).decs;
                Xr1 = Population(randi(end,1,N)).decs;
                Xr3 = Population(randi(end,1,N)).decs;
                P   = [Population,Archive];
                Xr2 = P(randi(end,1,N)).decs;
                CR  = randn(N,1).*sqrt(0.1) + MCR(randi(end,N,1));
                CR  = sort(CR);
                CR  = repmat(max(0,min(1,CR)),1,Problem.D);
                F   = min(1,trnd(1,N,1).*sqrt(0.1) + MF(randi(end,N,1)));
                while any(F<=0)
                    F(F<=0) = min(1,trnd(1,sum(F<=0),1).*sqrt(0.1) + MF(randi(end,sum(F<=0),1)));
                end
                F  = repmat(F,1,Problem.D);
                OP = arrayfun(@(S)find(rand<=cumsum(MOP),1),1:N);
                OP = arrayfun(@(S)find(OP==S),1:length(MOP),'UniformOutput',false);
                % Generate offspring
                PopDec = Population.decs;
                OffDec = PopDec;
                OffDec(OP{1},:) = PopDec(OP{1},:) + F(OP{1},:).*(Xp1(OP{1},:)-PopDec(OP{1},:)+Xr1(OP{1},:)-Xr2(OP{1},:));
                OffDec(OP{2},:) = PopDec(OP{2},:) + F(OP{2},:).*(Xp1(OP{2},:)-PopDec(OP{2},:)+Xr1(OP{2},:)-Xr3(OP{2},:));
                OffDec(OP{3},:) = F(OP{3},:).*(Xr1(OP{3},:)+Xp2(OP{3},:)-Xr3(OP{3},:));
                if rand < 0.4
                    Site = rand(size(CR)) > CR;
                    OffDec(Site) = PopDec(Site);
                else
                    p1 = randi(Problem.D,N,1);
                    p2 = arrayfun(@(S)find([rand(1,Problem.D),2]>CR(S,1),1),1:N);
                    for i = 1 : N
                        Site = [1:p1(i)-1,p1(i)+p2(i):Problem.D];
                        OffDec(i,Site) = PopDec(i,Site);
                    end
                end
                Offspring = Problem.Evaluation(OffDec);
                % Update the population and archive
                delta   = FitnessSingle(Population) - FitnessSingle(Offspring);
                replace = delta > 0;
                Archive = [Archive,Population(replace)];
                Archive = Archive(randperm(end,min(end,ceil(aRate*N))));
                Population(replace) = Offspring(replace);
                % Update CR, F, and probabilities of operators
                if any(replace)
                    w      = delta(replace)./sum(delta(replace));
                    MCR(k) = (w'*CR(replace,1).^2)./(w'*CR(replace,1));
                    MF(k)  = (w'*F(replace,1).^2)./(w'*F(replace,1));
                    k      = mod(k,length(MCR)) + 1;
                else
                    MCR(k) = 0.5;
                    MF(k)  = 0.5;
                end
                delta = max(0,delta./abs(FitnessSingle(Population)));
                if any(cellfun(@isempty,OP))
                	MOP = ones(1,3)/3;
                else
                    MOP = cellfun(@(S)mean(delta(S)),OP);
                    MOP = max(0.1,min(0.9,MOP./sum(MOP)));
                end
            end
        end
    end
end