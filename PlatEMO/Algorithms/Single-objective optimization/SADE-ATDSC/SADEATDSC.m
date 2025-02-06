classdef SADEATDSC < ALGORITHM
% <2022> <single> <real/integer> <expensive>
% Surrogate-assisted differential evolution with adaptation of training data selection criterion
% N         --- 100 --- Population size of DE
% F         --- 0.5 --- Scaling factor of DE
% CR        --- 0.9 --- Crossover rate of DE
% mut       ---   3 --- Mutation strategy ID 
% xov       ---   1 --- Crossover strategy ID
% delta     --- 0.2 --- Rate of validation data
% RBFN_lib  ---   2 --- RBFN library {newrbe, RBF + RBF_eval, rbfcreate + rbfinterp}
% kernel_id ---   2 --- RBFN kernel {'linear' 'cubic' 'gaussian' 'thinplate' 'multiquadric' 'inversemultiquadric'}

%------------------------------- Reference --------------------------------
% K. Nishihara and M. Nakata. Surrogate-assisted differential evolution
% with adaptation of training data selection criterion. Proceedings of the
% IEEE Symposium Series on Computational Intelligence, 2022, 1675–1682.
%--------------------------------------------------------------------------
   
% This function is generated to match with NKTLab experiment environment by
% Kei Nishihara. 2023/9/13
% The original code was written in Python3.

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [N,F,CR,mut,xov,delta,RBFN_lib,kernel_id,adapt_datacrtr,datacrtr_id,spread]...
                = Algorithm.ParameterSet(100,0.5,0.9,3,1,0.2,2,2,true,1,NaN);

            predefined = [1, 2, 3, 4];  % {1, 2, 3, 4, 5} : {All, Current Pop or Top, Recent, Neighbor, Near Best}
            
            switch RBFN_lib
                case 2
                    kernels = {'LN' 'CB' 'GA' 'TPS' 'MQ' 'IMQ'}; kernel = kernels{kernel_id};
                case 3
                    kernels = {'linear' 'cubic' 'gaussian' 'thinplate' 'multiquadric'}; kernel = kernels{kernel_id};
                otherwise
                    kernel = 'GA';
            end

            datacrtr_len = length(predefined);    
            D = Problem.D;

            %% Generate random population
            PopDec     = UniformPoint(N, D, 'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper - Problem.lower, N, 1) .* PopDec + repmat(Problem.lower, N, 1), [repmat(Algorithm.metric.runtime, N, 1)]);
            Arc        = Population;

            % Archive
            db_x = Population.decs; db_y = Population.objs;
            gen  = 1;

            %% Optimization
            while Algorithm.NotTerminated(Arc)

                % Get the best
                [B, I] = sort(db_y); 
                pop = db_x(I(1:N), :); fit = B(1:N);
                % [bf, bid] = min(fit); bx = pop(bid,:);
                cand = DEGetTrialVector(pop, fit, Problem.lower, Problem.upper, F, CR, mut, xov);
                
                % Build surrogate
                mtrcs = zeros(1, datacrtr_len); surrs = cell(1, datacrtr_len);
                if adapt_datacrtr
                    for j = 1 : datacrtr_len
                        [train_x, train_y] = GetTrain_SADEATDSC(predefined(j), db_x, db_y, pop, N, N);
                        
                        % Validate
                        c    = cvpartition(size(train_x, 1), 'HoldOut', delta);
                        idx  = c.test; train_Xs = train_x(~idx, :); test_Xs = train_x(idx, :); train_Fs = train_y(~idx, :); test_Fs = train_y(idx, :);
                        surr = GetSurrogate_SADEATDSC(train_Xs, train_Fs, RBFN_lib, spread, kernel);
                        surr = @(x) Predict_SADEATDSC(x, surr, RBFN_lib, kernel);
                        pred_Fs  = surr(test_Xs);
                        mtrcs(j) = sqrt(mean((test_Fs - pred_Fs).^2));  % RMSE
                        surrs{j} = surr;
                    end
                    [~, flag] = min(mtrcs); surr = surrs{flag};
                else
                    [train_x, train_y] = GetTrain_SADEATDSC(datacrtr_id, db_x, db_y, pop, N, N);
                    c    = cvpartition(size(train_x, 1), 'HoldOut', delta);
                    idx  = c.test; train_Xs = train_x(~idx, :); test_Xs = train_x(idx, :); train_Fs = train_y(~idx, :); test_Fs = train_y(idx, :);
                    surr = GetSurrogate_SADEATDSC(train_Xs, train_Fs, RBFN_lib, spread, kernel);
                    surr = @(x) Predict_SADEATDSC(x, surr, RBFN_lib, kernel);
                end

                % Screen
                cand_Fs      = surr(cand);
                [~, idx]     = min(cand_Fs);
                offspringDec = cand(idx, :);

                % evaluate solutions
                offspring = Problem.Evaluation(offspringDec, [Algorithm.metric.runtime]);
                Arc       = [Arc, offspring];
                Algorithm.NotTerminated(Arc);

                db_x = [db_x; offspring.dec];  db_y = [db_y; offspring.obj];
                [db_x, ia, ~] = unique(db_x, 'stable', 'rows'); db_y = db_y(ia);

                gen = gen + 1;
            end
        end
    end
end