classdef TriP < ALGORITHM
% <multi/many> <real/integer> <constrained>
% Tri-population based coevolutionary algorithm
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, L. Wang, and C. Lu, A tri-population based
% co-evolutionary framework for constrained multi-objective optimization
% problems, Swarm and Evolutionary Computation, 2022, 70: 101055.
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
            [delta,nr] = Algorithm.ParameterSet(0.9,2);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population  = Problem.Initialization();
            Z           = min(Population.objs,[],1);
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            alpha_c     = 2./(1+exp(1).^(-Problem.FE*10/(3*Problem.maxFE)))-1;
            para        = ceil(Problem.maxFE/Problem.N)/2 - ceil(Problem.FE/(3*Problem.N));

            %% Evaluate the Population
            Tc               = 0.9 * ceil(Problem.maxFE/Problem.N);
            last_gen         = 20;
            change_threshold = 1e-1;
            search_stage     = 1; % 1 for push stage,otherwise,it is in pull stage.
            max_change       = 1;
            epsilon_k        = 0;
            epsilon_0        = 0;
            cp               = 2;
            alpha            = 0.95;
            tao              = 0.05;
            ideal_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            nadir_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            arch             = archive(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen        = ceil(Problem.FE/(2*Problem.N));
                pop_cons   = Population.cons;
                cv         = overall_cv(pop_cons);
                population = [Population.decs,Population.objs,cv];
                rf         = sum(cv <= 1e-6) / Problem.N;
                ideal_points(gen,:) = Z;
                nadir_points(gen,:) = max(population(:,Problem.D + 1 : Problem.D + Problem.M),[],1);

                % The maximumrate of change of ideal and nadir points rk is calculated.
                if gen >= last_gen
                    max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen);
                end

                % The value of e(k) and the search strategy are set.
                if gen < Tc
                    if max_change <= change_threshold && search_stage == 1
                        search_stage = -1;
                        epsilon_0 = max(population(:,end),[],1);
                        epsilon_k = epsilon_0;
                    end
                    if search_stage == -1
                        epsilon_k =  update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp);
                    end
                else
                    epsilon_k = 0;
                end

                % For each solution
                for i = 1 : Problem.N
                    % Choose the parents
                    if rand < delta
                        P = B(i,randperm(size(B,2)));
                    else
                        P = randperm(Problem.N);
                    end

                    % Generate an offspring
                    Offspring = OperatorDE(Problem,Population(i),Population(P(1)),Population(P(2)));

                    % Update the ideal point
                    Z = min(Z,Offspring.obj);

                    g_old  = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                    g_new  = max(repmat(abs(Offspring.obj-Z),length(P),1).*W(P,:),[],2);
                    cv_old = overall_cv(Population(P).cons);
                    cv_new = overall_cv(Offspring.con) * ones(length(P),1);

                    if search_stage == 1 % Push Stage
                        Population(P(find(g_old>=g_new,nr))) = Offspring;
                    else  % Pull Stage  &&  An improved epsilon constraint-handling is employed to deal with constraints
                        Population(P(find(((g_old >= g_new) & (((cv_old <= epsilon_k) & (cv_new <= epsilon_k)) | (cv_old == cv_new)) | (cv_new < cv_old) ), nr))) = Offspring;
                    end
                end

                % Re-rank
                Population1 = Population1(randperm(Problem.N));
                Population2 = Population2(randperm(Problem.N));
                Lia   = ismember(Population2.objs,Population1.objs, 'rows');
                gamma = 1-sum(Lia)/Problem.N;            
                [Population1,FrontNo1] = EnvironmentalSelection(Population1,Problem.N,alpha_c);
                [Population2,FrontNo2] = EnvironmentalSelection_noCon(Population2,Problem.N,alpha_c,gamma,para);

                % Offspring Reproduction
                Population_all   = [Population1,Population2];
                RankSolution_all = [FrontNo1,FrontNo2];
                MatingPool = TournamentSelection(2,2*Problem.N,RankSolution_all);
                Offspring  = OperatorGAhalf(Problem,Population_all(MatingPool));

                % Environmental Selection
                alpha_c = 2./(1+exp(1).^(-Problem.FE*5/Problem.maxFE)) - 1; 
                para  = ceil(Problem.maxFE/Problem.N)/2 - ceil(Problem.FE/(2*Problem.N));
                [Population1,~] = EnvironmentalSelection([Population1,Offspring],Problem.N,alpha_c);
                [Population2,~] = EnvironmentalSelection_noCon([Population2,Offspring],Problem.N,alpha_c,gamma,para);

                % Output the non-dominated and feasible solutions.
                arch = archive([arch,Population],Problem.N);
                if Problem.FE >= Problem.maxFE
                    arch = archive([arch,Population1],Problem.N);
                    Population = arch;
                end
            end
        end
    end
end

function result = overall_cv(cv)
% The Overall Constraint Violation

    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end

function max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen)
% Calculate the Maximum Rate of Change

    delta_value = 1e-6 * ones(1,size(ideal_points,2));
    rz = abs((ideal_points(gen,:) - ideal_points(gen - last_gen + 1,:)) ./ max(ideal_points(gen - last_gen + 1,:),delta_value));
    nrz = abs((nadir_points(gen,:) - nadir_points(gen - last_gen + 1,:)) ./ max(nadir_points(gen - last_gen + 1,:),delta_value));
    max_change = max([rz, nrz]);
end

function result = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp)
    if rf < alpha
        result = (1 - tao) * epsilon_k;
    else
        result = epsilon_0 * ((1 - (gen / Tc)) ^ cp);
    end
end