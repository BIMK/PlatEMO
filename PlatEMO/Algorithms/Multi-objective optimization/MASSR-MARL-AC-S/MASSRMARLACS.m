classdef MASSRMARLACS < ALGORITHM
% <2026> <multi> <real/binary> <large/none> <sparse>
% Sparse multi-objective optimization with MASSR states

%------------------------------- Reference --------------------------------
% X. Liu, Y. Tian, S. Yang, Z. Sheng, and X. Zhang. A multi-agent
% self-supervised state representation framework for automated algorithm
% configuration. IEEE Transactions on Evolutionary Computation, 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Estimate variable importance and initialize search
            N                                        = Problem.N;
            [Population,Dec,Mask,Fitness,Importance] = MASSRMARLACSInitialize(Problem);

            %% Initialize representation learning and sparse statistics
            warmFE    = ceil(0.2*Problem.maxFE);
            reference = MASSRMARLACSReference(Problem);
            encoder   = MASSREncoder(Problem.D+Problem.M,'SPARSE');
            encoder.collect(MASSRFeatures.observation(Population,Problem), ...
            MASSRFeatures.guidance(Population,[],Problem,0,0,0));
            controller = [];
            sv         = 0.5*ones(1,Problem.D);
            pv         = 0.5*ones(1,Problem.D);
            fv         = zeros(1,Problem.D);
            Memory     = [];
            totalElite = 0;
            stagnation = 0;
            window     = max(1,0.1*Problem.maxFE/N);
            score      = HV(Population,reference);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Train MASSR once at the random-collection boundary
                if ~encoder.Ready && Problem.FE == warmFE
                    encoder.fit();
                    controller = MASSRMARLACSPolicy();
                end

                % Select guiding, ratio, and variation actions
                if encoder.Ready
                    state           = encoder.encode(MASSRFeatures.observation(Population,Problem));
                    [action,vector] = controller.act(state,Problem.FE/Problem.maxFE);
                    limit           = Problem.maxFE;
                else
                    action = [randi(3);rand;randi(2)];
                    limit  = warmFE;
                end

                % Construct the selected variable-guiding vector
                LastPopulation = Population;
                switch action(1)
                    case 1
                        weights = Importance;
                    case 2
                        weights = sum(~Mask,1);
                    case 3
                        weights = rand(1,Problem.D);
                end

                % Update elite activation frequencies and variable spread
                progress = Problem.FE/Problem.maxFE;
                elite    = Fitness < 1;
                if any(elite)
                    number     = sum(elite);
                    sv         = (totalElite*sv+sum(Mask(elite,:),1))/(totalElite+number);
                    totalElite = totalElite + number;
                    if progress < 0.618
                        fv = std(Dec(elite,:),0,1);
                    end
                end
                if progress < 0.618
                    pv = pv.*(1-sv)*sqrt(progress) + pv;
                end

                % Generate masked offspring within the current FE limit
                count            = min(N,limit-Problem.FE);
                parents          = TournamentSelection(2,2*count,Fitness);
                [OffDec,OffMask] = MASSRMARLACSOperator(Problem,Dec(parents,:),Mask(parents,:), ...
                    weights,Mask,action(3),Memory,pv,fv,progress);
                Offspring        = Problem.Evaluation(OffDec.*OffMask);

                % Select the population size from the normalized ratio action
                sizeNext                      = ceil(N*(1.5+0.5*action(2)));
                candidates                    = [Population,Offspring];
                [Population,Dec,Mask,Fitness] = MASSRMARLACSSelection(candidates, ...
                    [Dec;OffDec],[Mask;OffMask],sizeNext);

                % Compute the unscaled HV change as the team reward
                previous    = score;
                score       = HV(Population,reference);
                reward      = score - previous;
                improvement = reward/max(abs(previous),1e-12);
                if improvement > 1e-6
                    stagnation = 0;
                else
                    stagnation = stagnation + 1;
                end
                replacement = mean(~ismember(Population.decs,LastPopulation.decs,'rows'));

                % Update guiding history and store the learning sample
                Memory = [Memory,action(1)];
                Memory = Memory(max(1,end-49):end);
                X      = MASSRFeatures.observation(Population,Problem);
                if encoder.Ready
                    controller.remember(state,vector,reward,encoder.encode(X),Problem.FE>=Problem.maxFE);
                    controller.train();
                else
                    encoder.collect(X,MASSRFeatures.guidance(Population,LastPopulation,Problem, ...
                    improvement,stagnation/window,replacement));
                end
            end
        end
    end
end