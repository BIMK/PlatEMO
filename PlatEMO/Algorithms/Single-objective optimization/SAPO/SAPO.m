classdef SAPO < ALGORITHM
% <2024> <single> <real/integer> <constrained> <expensive>
% Surrogate-assisted partial optimization
% N             --- 100   --- Population size of DE
% F             --- 0.5   --- Scaling factor of DE
% CR            --- 0.9   --- Crossover rate of DE
% mut1          --- 1     --- Mutation strategy ID 
% mut2          --- 3     --- Mutation strategy ID 
% xov           --- 1     --- Crossover strategy ID
% initsize1     --- 100   --- Initial sampling size for LHS
% initsize2     --- 200   --- Initial sampling size for LHS
% change_data_n --- 2     --- {0, 1, 2} : whether change datasize {0: False, 1:True-fix, 2:True-nD}
% ch_d_n_thres  --- 100   --- Threshold to change datasize if change_data_n = 1
% data_n        --- 100   --- The number of neighbor if change_data_n = 0
% data_n_2      --- 200   --- The number of data size if change_data_n = 1
% data_times    --- 5     --- Parameter to define data size if change_data_n = 2

%------------------------------- Reference --------------------------------
% K. Nishihara and M. Nakata. A surrogate-assisted partial optimization for
% expensive constrained optimization problems. Proceedings of the
% International Conference on Parallel Problem Solving from Nature, 2024,
% 391–407.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [N,F,CR,mut1,mut2,xov,initsize1,initsize2,change_data_n,ch_d_n_thres,data_n,data_n_2,data_times,kernel]...
                = Algorithm.ParameterSet(100,0.5,0.9,1,3,1,100,200,2,100,100,200,5,'cubic');

            D = Problem.D;

            if D < 100
                init_size = initsize1;
            else
                init_size = initsize2;
            end

            %% Generate random population
            PopDec     = UniformPoint(init_size, D, 'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper - Problem.lower, init_size, 1) .* PopDec + repmat(Problem.lower, init_size, 1), [repmat(Algorithm.metric.runtime, init_size, 1)]);
            CV         = sum(max(0,Population.cons),2); 
            for i = 1 : size(CV,1)
                Population(i).add = [Population(i).add, CV(i)];
            end
            Arc = Population;

            % Archive
            db_x = Population.decs; db_y = Population.objs; db_c = Population.cons; 
            db_v = max(0, db_c); db_cv = sum(db_v, 2); db_feasible = (db_cv == 0);

            ncon = size(Population.cons, 2);
            gen  = 1;

            % Note that both equality and inequality constraints are denoted as g

            %% Optimization
            while Algorithm.NotTerminated(Arc)
                
                % Get the best
                [~, ~, ~, ~, bf] = Findbest_SAPO(db_x, db_y, db_c, db_cv);

                % Find indices of feasible but better than bf solutions
                ind_infea_better = find(and(~db_feasible, (db_y <= bf)));
                len = size(db_x, 1);
                ind_worse = setdiff((1:len).', ind_infea_better);   % complement

                switch rem(gen, 2)
                    
                    case 1  % f, gm -> gn
                        for n = 1 : ncon
                            
                            [pop, fit, con, train_x, train_y, train_c] = Select_fg4g_RtnTrain_SAPO(db_x, db_y, db_c, ind_infea_better, ind_worse, ncon, N, n, change_data_n, ch_d_n_thres, data_n, data_n_2, data_times);

                            [~, bid] = Findbest_f2g_SAPO(n, pop, fit, con);
                            
                            if size(pop, 1) < 4
                                continue
                            end

                            cand = [
                                    DEGetTrialVectorConstrained(pop, [], [], F, CR, mut1, xov, bid); ...
                                    DEGetTrialVectorConstrained(pop, [], [], F, CR, mut2, xov, bid) ...
                                    ];
                            cand = Repair_SAPO(cand, Problem.lower, Problem.upper);

                            % Construct surrogate models for each f and g
                            para = RBFCreate_SAPO(train_x, [train_y, train_c], kernel);
                            
                            cand_return = RBFInterp_SAPO(cand, para);
                            cand_fit    = cand_return(:, 1); cand_con = cand_return(:, 2:end);
                            cand_bx     = Findbest_f2g_SAPO(n, cand, cand_fit, cand_con);
        
                            % evaluate solution
                            Offspring = Problem.Evaluation(cand_bx, [Algorithm.metric.runtime]);
                            CV        = sum(max(0,Offspring.con), 2); 
                            Offspring.add = [Offspring.add, CV];
                            Arc       = [Arc, Offspring];
                            Algorithm.NotTerminated(Arc);
        
                            db_x = [db_x; Offspring.dec]; db_y = [db_y; Offspring.obj]; db_c = [db_c; Offspring.con]; 
                            v    = max(0, Offspring.con); db_v = [db_v; v]; 
                            cv   = sum(v, 2); db_cv = [db_cv; cv];
                            db_feasible = [db_feasible; (cv == 0)];
                        end


                    case 0  % gm -> f

                        for n = 1 : ncon
                            [pop, fit, cv, train_x, train_y, train_c] = Select_g4f_half_RtnTrain_SAPO(db_x, db_y, db_c, ind_infea_better, ind_worse, ncon, N, db_feasible, len, change_data_n, ch_d_n_thres, data_n, data_n_2, data_times);
                            
                            [~, bid] = Findbest_g2f_SAPO(pop, fit, cv);  % Same to the Feasibility Rule
                            
                            if size(pop, 1) < 4
                                continue
                            end
                            
                            cand = [
                                    DEGetTrialVectorConstrained(pop, [], [], F, CR, mut1, xov, bid); ...
                                    DEGetTrialVectorConstrained(pop, [], [], F, CR, mut2, xov, bid) ...
                                    ];
                            cand = Repair_SAPO(cand, Problem.lower, Problem.upper);
    
                            % Construct surrogate models for each f and g
                            para = RBFCreate_SAPO(train_x, [train_y, train_c], kernel);
                                    
                            cand_return = RBFInterp_SAPO(cand, para);
                            cand_fit    = cand_return(:, 1); cand_con = cand_return(:, 2:end);
                            cand_v      = max(0, cand_con); cand_cv = sum(cand_v, 2);
                            cand_bx     = Findbest_g2f_SAPO(cand, cand_fit, cand_cv);  % Same to the Feasibility Rule
        
                            % evaluate solution
                            Offspring = Problem.Evaluation(cand_bx, [Algorithm.metric.runtime]);
                            CV        = sum(max(0,Offspring.con), 2); 
                            Offspring.add = [Offspring.add, CV];
                            Arc       = [Arc, Offspring];
                            Algorithm.NotTerminated(Arc);
        
                            db_x = [db_x; Offspring.dec]; db_y = [db_y; Offspring.obj]; db_c = [db_c; Offspring.con]; 
                            v    = max(0, Offspring.con); db_v = [db_v; v]; 
                            cv   = sum(v, 2); db_cv = [db_cv; cv];
                            db_feasible = [db_feasible; (cv == 0)];
                        end
                end
                gen = gen + 1;        
            end
        end
    end
end