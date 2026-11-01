classdef MASSRCMODRL < ALGORITHM
% <2026> <multi> <real> <constrained>
% Constrained multi-objective optimization with MASSR states

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
            %% Generate the two populations and archive
            N                      = Problem.N;
            N2                     = ceil(N/2);
            Population1            = Problem.Initialization(N);
            Population2            = Problem.Initialization(N2);
            [Archive,~,~]          = MASSRCMODRLSelection([Population1,Population2],N,1,0);
            [Population1,Fitness1] = MASSRCMODRLSelection(Population1,N,1,0);
            [Population2,Fitness2] = MASSRCMODRLSelection(Population2,N2,0,0);

            %% Initialize representation learning
            warmFE     = ceil(0.2*Problem.maxFE);
            encoder    = MASSREncoder(Problem.D+Problem.M,'CMODRL');
            controller = [];
            stagnation = 0;
            window     = max(1,0.1*Problem.maxFE/(N+N2));
            encoder.collect(MASSRFeatures.observation(Population1,Problem), ...
            MASSRFeatures.guidance(Population1,[],Problem,0,0,0));

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                % Train MASSR once at the random-collection boundary
                if ~encoder.Ready && Problem.FE == warmFE
                    encoder.fit();
                    controller = MASSRCMODRLPolicy();
                end

                % Select random or learned configuration actions
                if encoder.Ready
                    state  = encoder.encode(MASSRFeatures.observation(Population1,Problem));
                    action = controller.act(state);
                    limit  = Problem.maxFE;
                else
                    action = [rand;randi(2)];
                    limit  = warmFE;
                end

                % Generate offspring within the current FE limit
                LastPopulation = Population1;
                LastArchive    = Archive;
                averageCV      = mean(sum(max(0,Population1.cons),2));
                epsilon        = action(1)*averageCV;
                count          = min(N+N2,limit-Problem.FE);
                n1             = min(N,count);
                n2             = count - n1;
                Offspring1     = MASSRCMODRLOffspring(Problem,Population1,Fitness1,action(2),n1);
                Offspring2     = MASSRCMODRLOffspring(Problem,Population2,Fitness2,action(2),n2);
                children       = [Offspring1,Offspring2];

                % Update the archive and both populations
                [Archive,~,~]                   = MASSRCMODRLSelection([Archive,children],N,1,0);
                [Population1,Fitness1,selected] = MASSRCMODRLSelection([Population1,children],N,2,epsilon);
                [Population2,Fitness2]          = MASSRCMODRLSelection([Population2,children],N2,0,0);

                % Compute the shared archive-improvement reward
                old         = LastArchive.objs;
                now         = Archive.objs;
                convergence = mean(sum(old,2)) - mean(sum(now,2));
                diversity   = sum(max(now,[],1)-min(now,[],1)) - sum(max(old,[],1)-min(old,[],1));
                feasibility = mean(sum(max(0,LastArchive.cons),2)) - mean(sum(max(0,Archive.cons),2));
                reward      = convergence + diversity + feasibility;
                improvement = MASSRFeatures.change(mean(sum(old,2)),mean(sum(now,2)));
                if improvement > 1e-6
                    stagnation = 0;
                else
                    stagnation = stagnation + 1;
                end

                % Store a policy transition or a representation sample
                replacement = sum(selected(N+1:end))/N;
                X           = MASSRFeatures.observation(Population1,Problem);
                if encoder.Ready
                    controller.remember(state,action,reward,encoder.encode(X),Problem.FE>=Problem.maxFE);
                    controller.train();
                else
                    encoder.collect(X,MASSRFeatures.guidance(Population1,LastPopulation, ...
                    Problem,improvement,stagnation/window,replacement));
                end
            end
        end
    end
end