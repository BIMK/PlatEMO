classdef PPS < ALGORITHM
% <multi/many> <real/integer> <constrained>
% Push and pull search algorithm
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Li, C. Wei, Q. Zhang, K. Deb, and E. Goodman,
% Push and pull search for solving constrained multi-objective optimization
% problems, Swarm and Evolutionary Computation, 2019, 44(2): 665-679.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenji Li

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
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);

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
                gen        = ceil(Problem.FE/Problem.N);
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
                        epsilon_k = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp);
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

                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),length(P),1).*W(P,:),[],2);
                    cv_old = overall_cv(Population(P).cons);
                    cv_new = overall_cv(Offspring.con) * ones(length(P),1);

                    if search_stage == 1 % Push Stage
                        Population(P(find(g_old>=g_new,nr))) = Offspring;
                    else  % Pull Stage  &&  An improved epsilon constraint-handling is employed to deal with constraints
                        Population(P(find(((g_old >= g_new) & (((cv_old <= epsilon_k) & (cv_new <= epsilon_k)) | (cv_old == cv_new)) | (cv_new < cv_old) ), nr))) = Offspring;
                    end
                end

                % Output the non-dominated and feasible solutions.
                arch = archive([arch,Population],Problem.N);
                if Problem.FE >= Problem.maxFE
                    Population = arch;
                end
            end
        end
    end
end

% The Overall Constraint Violation
function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end

% Calculate the Maximum Rate of Change
function max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen)
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