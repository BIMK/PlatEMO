classdef CMOEMT < ALGORITHM
% <multi> <real> <constrained>
% Constrained multi-objective optimization based on evolutionary multitasking optimization
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, L. Wang, and L. Gao, Constrained multi-objective
% optimization via multitasking and knowledge transfer, IEEE Transactions
% on Evolutionary Computation, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

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
            Population{1} = Problem.Initialization();
            Population{2} = Problem.Initialization();
            Z             = min(Population{2}.objs,[],1);
            Population{3} = Problem.Initialization();
            LengthO{1} = Problem.N/2; LengthO{2} = Problem.N/2; LengthO{3} = Problem.N/2;
            Fitness{1}    = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{3}    = CalFitness(Population{3}.objs);
            
            %% Evaluate the Population
            Tc               = 0.9 * ceil(Problem.maxFE/Problem.N);
            last_gen         = 20;
            change_threshold = 1e-1;
            search_stage     = 1; % 1 for push stage,otherwise,it is in pull stage.
            max_change       = 1;
            epsilon_k        = 0;
            epsilon_0        = 0;
            cp               = 2;
            alpha1            = 0.95;
            tao              = 0.05;
            ideal_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            nadir_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            
            %% Optimization
            while Algorithm.NotTerminated(Population{1})
                gen        = ceil(Problem.FE/(2*Problem.N));
                pop_cons   = Population{2}.cons;
                cv         = overall_cv(pop_cons);
                population = [Population{2}.decs,Population{2}.objs,cv];
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
                        epsilon_k =  update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha1,gen,Tc,cp);
                    end
                else
                    epsilon_k = 0;
                end
                
                if Problem.FE < Problem.maxFE/2
                    %non transfer
                    
                    % Offspring Reproduction
                    MatingPool{1} = TournamentSelection(2,2*LengthO{1},Fitness{1});
                    Offspring{1}  = OperatorGAhalf(Problem,Population{1}(MatingPool{1}));
                    Offspring{2} = [];
                    
                    for subgeneration = 1 : 5
                        Bounday = find(sum(W<1e-3,2)==Problem.M-1)';
                        Bounday = [Bounday,floor(length(W)/2)];
                        I = [Bounday,randi(length(W),1,floor(Problem.N/5)-length(Bounday))];
                        for j = 1 : length(I)
                            i = I(j);
                            if rand < delta
                                P = B(i,randperm(size(B,2)));
                            else
                                P = randperm(Problem.N);
                            end
                            
                            % Generate an offspring
                            offspring = OperatorDE(Problem,Population{2}(i),Population{2}(P(1)),Population{2}(P(2)));
                            Offspring{2} = [Offspring{2},offspring];
                            
                            % Update the ideal point
                            Z = min(Z,offspring.obj);
                            
                            % TCH approach
                            g_old = max(abs(Population{2}(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                            g_new = max(repmat(abs(offspring.obj-Z),length(P),1).*W(P,:),[],2);
                            cv_old = overall_cv(Population{2}(P).cons);
                            cv_new = overall_cv(offspring.con) * ones(length(P),1);
                            
                            if search_stage == 1 % Push Stage
                                Population{2}(P(find(g_old>=g_new,nr))) = offspring;
                            else  % Pull Stage  &&  An improved epsilon constraint-handling is employed to deal with constraints
                                Population{2}(P(find(((g_old >= g_new) & (((cv_old <= epsilon_k) & (cv_new <= epsilon_k)) | (cv_old == cv_new)) | (cv_new < cv_old) ), nr))) = offspring;
                            end
                            
                        end
                    end
                    
                    MatingPool{3} = TournamentSelection(2,2*LengthO{3},Fitness{3});
                    Offspring{3}  = OperatorGAhalf(Problem,Population{3}(MatingPool{3}));
                    
                    % Environmental Selection
                    [Population{1},Fitness{1},~] = EnvironmentalSelectionT1([Population{1},Offspring{1:3}],Problem.N);     
                    [Population{3},Fitness{3},~] = EnvironmentalSelectionT3([Population{3},Offspring{1:3}],Problem.N);
                    
                else
                    %transfer
                    % Offspring Reproduction
                    MatingPool{1} = TournamentSelection(2,2*LengthO{1},Fitness{1});
                    Offspring{1}  = OperatorGAhalf(Problem,Population{1}(MatingPool{1}));
                    Offspring{2}  = [];

                    for subgeneration = 1 : 5
                        Bounday = find(sum(W<1e-3,2)==Problem.M-1)';
                        Bounday = [Bounday,floor(length(W)/2)];
                        I = [Bounday,randi(length(W),1,floor(Problem.N/5)-length(Bounday))];
                        for j = 1 : length(I)
                            i = I(j);
                            
                            %                    for i = 1 : Problem.N
                            % Choose the parents
                            if rand < delta
                                P = B(i,randperm(size(B,2)));
                            else
                                P = randperm(Problem.N);
                            end
                            
                            % Generate an offspring
                            offspring = OperatorDE(Problem,Population{2}(i),Population{2}(P(1)),Population{2}(P(2)));
                            Offspring{2} = [Offspring{2},offspring];
                            
                            % Update the ideal point
                            Z = min(Z,offspring.obj);
                            
                            % TCH approach
                            g_old = max(abs(Population{2}(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                            g_new = max(repmat(abs(offspring.obj-Z),length(P),1).*W(P,:),[],2);
                            cv_old = overall_cv(Population{2}(P).cons);
                            cv_new = overall_cv(offspring.con) * ones(length(P),1);
                        
                            if search_stage == 1 % Push Stage
                                Population{2}(P(find(g_old>=g_new,nr))) = offspring;
                            else  % Pull Stage  &&  An improved epsilon constraint-handling is employed to deal with constraints
                                Population{2}(P(find(((g_old >= g_new) & (((cv_old <= epsilon_k) & (cv_new <= epsilon_k)) | (cv_old == cv_new)) | (cv_new < cv_old) ), nr))) = offspring;
                            end
                            
                        end
                    end
                    
                    MatingPool{3} = TournamentSelection(2,2*LengthO{3},Fitness{3});
                    Offspring{3}  = OperatorGAhalf(Problem,Population{3}(MatingPool{3}));
                    
                    % Environmental Selection
                    [Population{1},Fitness{1},~] = EnvironmentalSelectionT1([Population{1:3},Offspring{1:3}],Problem.N);     
                    [Population{3},Fitness{3},~] = EnvironmentalSelectionT3([Population{1:3},Offspring{1:3}],Problem.N);
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