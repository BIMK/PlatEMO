classdef CSEMT < ALGORITHM
% <2024> <multi> <real/integer/label/binary/permutation> <constrained>
% Constraints separation based evolutionary multitasking

%------------------------------- Reference --------------------------------
% K. Qiao, J. Liang, K. Yu, X. Ban, C. Yue, B. Qu, and P. N. Suganthan.
% Constraints separation based evolutionary multitasking for constrained
% multi-objective optimization problems. IEEE/CAA Journal of Automatica
% Sinica, 2024, 11(8): 1819-1835.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    methods
        function main(Algorithm,Problem)
            %% Generate random population 
            Population{1} = Problem.Initialization();   % main task considering all constraints
            number_cons   = size(Population{1}.cons,2); % the number of constraints
            Fitness{1}    = ICalFitness(Population{1}.objs,Population{1}.cons,1:number_cons);

            max_tasks = number_cons+1;

            for i = 1 : number_cons
                target{i}     = i; % the index of constraint
                archive{i}    = Problem.Initialization();
                ar_Fitness{i} = ICalFitness(archive{i}.objs,archive{i}.cons,[target{i}]);
            end
            target{max_tasks}     = [];
            archive{max_tasks}    = Problem.Initialization();
            ar_Fitness{max_tasks} = ICalFitness(archive{max_tasks}.objs,archive{max_tasks}.cons,[target{i}]);

            % Parameters for judging the stage switching
            gen              = 0;
            last_gen         = 100;
            change_threshold = 0.5;
            change_rate      = zeros(ceil(Problem.FE/Problem.N),Problem.M);

            cnt  = 0;   % evolutionary generation
            flag = 0;   % evolutionary stage

            is_flag    = [];
            is_flag_ar = [];
            sd         = [];
            beta       = 0.8;

            %% Optimization
            while Algorithm.NotTerminated(Population{1})
                cnt = cnt + 1;
                if flag == 0
                    change_rate = Normalization(archive{max_tasks},change_rate,cnt);
                    if Convertion(change_rate,cnt,gen,last_gen,change_threshold)
                        flag = 1;
                        for i = 1 : max_tasks-1
                            for j = i+1 : max_tasks
                                if i~= j
                                    if i == max_tasks
                                        c1=1;
                                    else
                                        c1 = 0;
                                    end
                                    if j == max_tasks
                                        c2 = 1;
                                    else
                                        c2 = 0;
                                    end

                                    is_flag(i,j) = Analysis_of_relationship(archive{i},archive{j},c1,c2,target{i},target{j},beta);

                                    if is_flag(i,j) == 2
                                        is_flag(j,i) =2;
                                    elseif is_flag(i,j) == 1
                                        is_flag(j,i) =0;
                                    elseif is_flag(i,j) ==0
                                        is_flag(j,i) =1;
                                    end
                                end
                            end
                        end
                        sd = [];
                        domi = cell(1,max_tasks);
                        for i = 1 : max_tasks
                            for j = 1 : max_tasks
                                if j~= i
                                    if is_flag(i,j)  == 0
                                        domi{i} = [domi{i},j];
                                    end
                                end
                            end
                            if isempty(domi{i})
                                sd = [sd,i];
                            end
                        end

                        for i = 1 : length(sd)
                            is_flag_ar(i,1) = Analysis_of_relationship2(archive{sd(i)},Population{1},target{sd(i)},[1:number_cons],beta);
                        end

                        % retained auxiliary populations
                        new_archive = archive(sd);
                        new_archive_fitness = ar_Fitness(sd);
                        target      = target(sd);
                        max_tasks   = length(sd);
                        archive     = new_archive;
                        ar_Fitness  = new_archive_fitness;
                    end
                end

                MatingPool1  = TournamentSelection(2,Problem.N,Fitness{1});
                Offspring{1} = OperatorGAhalf(Problem,Population{1}(MatingPool1));

                main_total_Off = Offspring{1};
                fu_total_Off   = cell(max_tasks,1);

                for i = 1 : max_tasks
                    MatingPool2 = TournamentSelection(2,Problem.N,ar_Fitness{i});
                    archive_Offspring{i} = OperatorGAhalf(Problem,archive{i}(MatingPool2));

                    if flag == 0 % the first stage
                        main_total_Off = [main_total_Off,archive_Offspring{i}];
                    else % the second stage
                        [~,~,Next] = EnvironmentalSelection( [archive{i},archive_Offspring{i}],Problem.N,[1:number_cons]);
                        succ_rate(i,1) = (sum(Next(1:length(archive{i})))/length(archive{i})) - (sum(Next(length(archive{i})+1:end))/length(archive_Offspring{i}));

                        if is_flag_ar(i,1) ~= 0
                            if succ_rate(i,1)>0
                                rand_number = randperm(Problem.N);
                                main_total_Off = [main_total_Off,archive{i}(rand_number(1:Problem.N/2))];
                            else
                                main_total_Off = [main_total_Off,archive_Offspring{i}];
                            end
                        else
                            rand_number    = randperm(Problem.N);
                            main_total_Off = [main_total_Off,archive{i}(rand_number(1:Problem.N/2))];
                            main_total_Off = [main_total_Off,archive_Offspring{i}];
                        end

                        [~,~,Next]     = EnvironmentalSelection( [Population{1},Offspring{1}],Problem.N,[target{i}]);
                        succ_rate(i,2) = (sum(Next(1:length(Population{1})))/length(Population{1})) - (sum(Next(length(Population{1})+1:end))/length(Offspring{1}));
                    end
                end

                for i = 1 : max_tasks
                    if flag == 0
                        fu_total_Off{i} = Offspring{1};
                    else
                        if is_flag_ar(i,1) ~= 0
                            if succ_rate(i,2) > 0
                                rand_number     = randperm(Problem.N);
                                fu_total_Off{i} = [fu_total_Off{i},Population{1}(rand_number(1:Problem.N/2))];
                            else
                                fu_total_Off{i} = [fu_total_Off{i},Offspring{1}];
                            end
                        else
                            rand_number     = randperm(Problem.N);
                            fu_total_Off{i} = [fu_total_Off{i},Population{1}(rand_number(1:Problem.N/2))];

                            fu_total_Off{i} = [fu_total_Off{i},Offspring{1}];
                        end
                    end
                    for j = 1 : max_tasks
                        fu_total_Off{i} = [fu_total_Off{i},archive_Offspring{j}];
                    end
                end

                % Environmental selection
                [Population{1},Fitness{1}] = EnvironmentalSelection([Population{1},main_total_Off],Problem.N,[1:number_cons]);
                for i = 1 : max_tasks
                    [archive{i}, ar_Fitness{i}] = EnvironmentalSelection([archive{i},fu_total_Off{i}],Problem.N,[target{i}]);
                end
            end
        end
    end
end