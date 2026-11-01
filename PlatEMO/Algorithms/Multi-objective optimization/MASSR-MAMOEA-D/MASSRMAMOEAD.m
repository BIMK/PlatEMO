classdef MASSRMAMOEAD < ALGORITHM
% <2026> <multi> <real>
% Multi-agent self-supervised state representation with MOEA/D

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

    properties(SetAccess = private)
        ReferenceFront;    % Fixed Pareto reference set used for IGD rewards
    end

    methods
        function main(Algorithm,Problem)
            %% Initialize decomposition and reference settings
            [W,Problem.N]            = UniformPoint(Problem.N,Problem.M);
            N                        = Problem.N;
            warmFE                   = ceil(0.2*Problem.maxFE);
            Algorithm.ReferenceFront = MASSRMAMOEADReference(Problem);

            %% Generate the initial population and archive
            Population      = Problem.Initialization();
            Archive         = Population;
            Z               = min(Population.objs,[],1);
            [~,neighbors]   = sort(MASSRFeatures.distance(W,W),2);
            interval        = max(1,floor(0.05*Problem.maxFE/N));
            cooling         = interval;
            Toptions        = [15,20,25,30];
            Foptions        = [0.4,0.5,0.6,0.7];
            score           = IGD(Population,Algorithm.ReferenceFront);
            [~,tracker]     = MASSRMAMOEADReward(score,[]);
            stagnation      = 0;
            stagnationScale = max(1,0.1*Problem.maxFE/N);

            %% Initialize representation learning
            encoder = MASSREncoder(Problem.D+Problem.M,'MAMOEAD');
            encoder.collect(MASSRFeatures.observation(Population,Problem), ...
            MASSRFeatures.guidance(Population,[],Problem,0,0,0));
            controller = [];

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Train MASSR once at the random-collection boundary
                if ~encoder.Ready && Problem.FE == warmFE
                    encoder.fit();
                    controller = MASSRMAMOEADPolicy();
                end

                % Select the four configuration actions
                available = availableActions(cooling,Problem.FE,Problem.maxFE);
                if encoder.Ready
                    state    = encoder.encode(MASSRFeatures.observation(Population,Problem));
                    progress = (Problem.FE-warmFE)/max(Problem.maxFE-warmFE-N,1);
                    action   = controller.act(state,available,progress);
                    limit    = Problem.maxFE;
                else
                    action = zeros(4,1);
                    for k = 1 : 4
                        valid     = find(available(:,k));
                        action(k) = valid(randi(numel(valid)));
                    end
                    limit = warmFE;
                end

                % Adjust weights when the selected action is available
                LastPopulation = Population;
                if action(1) == 2
                    [Population,W] = MASSRMAMOEADWeights(Population,W,Z,Archive);
                    [~,neighbors]  = sort(MASSRFeatures.distance(W,W),2);
                    cooling        = interval;
                end

                % Generate offspring and update decomposition subproblems
                T            = min(Toptions(action(2)),N);
                count        = min(N,limit-Problem.FE);
                order        = randperm(N);
                Offspring    = Population([]);
                replacements = 0;
                for i = order(1:count)
                    if rand < 0.8
                        mating = neighbors(i,1:T);
                    else
                        mating = 1 : N;
                    end
                    pool             = mating(mating~=i);
                    parent           = pool(randperm(numel(pool),5));
                    child            = MASSRMAMOEADOperator(Problem,Population(i).dec, ...
                        Population(parent).decs,action(3),Foptions(action(4)));
                    Offspring(end+1) = child;
                    Z                = min(Z,child.obj);
                    candidates       = mating(randperm(numel(mating)));
                    old              = max(abs(Population(candidates).objs-Z).*max(W(candidates,:),1e-6),[],2);
                    new              = max(abs(child.obj-Z).*max(W(candidates,:),1e-6),[],2);
                    selected         = candidates(find(new<old,2));
                    if ~isempty(selected)
                        Population(selected) = child;
                    end
                    replacements = replacements + numel(selected);
                end

                % Update the archive and best-so-far IGD reward
                Archive          = updateArchive(Archive,Offspring,floor(1.5*N));
                cooling          = max(cooling-1,0);
                previousScore    = score;
                score            = IGD(Population,Algorithm.ReferenceFront);
                [reward,tracker] = MASSRMAMOEADReward(score,tracker);
                improvement      = (previousScore-score)/max(abs(previousScore),1e-12);
                if improvement > 1e-6
                    stagnation = 0;
                else
                    stagnation = stagnation + 1;
                end

                % Store a policy transition or a representation sample
                Xnext = MASSRFeatures.observation(Population,Problem);
                if encoder.Ready
                    nextState     = encoder.encode(Xnext);
                    nextAvailable = availableActions(cooling,Problem.FE,Problem.maxFE);
                    controller.remember(state,action,reward,nextState, ...
                    Problem.FE>=Problem.maxFE,nextAvailable);
                    controller.train();
                else
                    encoder.collect(Xnext,MASSRFeatures.guidance(Population,LastPopulation, ...
                    Problem,improvement,stagnation/stagnationScale,min(replacements/N,1)));
                end
            end
        end
    end
end

function mask = availableActions(cooling,FE,maxFE)
%availableActions - Mask unavailable MOEA/D weight adjustments.
%
%   Weight adjustment is enabled only after cooldown and before 90% of the
%   evaluation budget. Other action components remain available.

    mask        = true(4,4);
    mask(3:4,1) = false;
    mask(2,1)   = cooling == 0 && FE < 0.9*maxFE;
end

function Archive = updateArchive(Archive,Offspring,capacity)
%updateArchive - Maintain a bounded nondominated archive.
%
%   Keep unique nondominated objective vectors. If capacity is exceeded,
%   repeatedly remove the solution with the smallest local distance product.

    Archive        = [Archive,Offspring];
    rank           = NDSort(Archive.objs,1);
    Archive        = Archive(rank==1);
    [~,uniqueRows] = unique(Archive.objs,'rows','stable');
    Archive        = Archive(uniqueRows);
    if numel(Archive) <= capacity
        return;
    end
    D                    = MASSRFeatures.distance(Archive.objs,Archive.objs);
    D(1:size(D,1)+1:end) = inf;
    while numel(Archive) > capacity
        ordered         = sort(D,2);
        sparsity        = prod(ordered(:,1:min(size(D,1)-1,size(Archive.objs,2))),2);
        [~,remove]      = min(sparsity);
        Archive(remove) = [];
        D(remove,:)     = [];
        D(:,remove)     = [];
    end
end