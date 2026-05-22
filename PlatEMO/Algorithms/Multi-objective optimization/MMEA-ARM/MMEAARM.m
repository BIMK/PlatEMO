classdef MMEAARM < ALGORITHM
% <2026> <multi> <real/integer/binary> <large/none> <multimodal> <sparse>
% Adaptive resource management based MMEA

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, and Y. Zhang. An adaptive resource management based
% evolutionary algorithm for large-scale multi-modal multi-objective
% optimization. Applied Soft Computing, 2026, 193: 114853.
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
            %% Population initialization
            Dec  = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1)); 
            Mask = rand(Problem.N, Problem.D) < 0.5;
            Population  = Problem.Evaluation(Dec.*Mask);
            K           = 10;  
            Masks       = cell(1,K);
            Decs        = cell(1,K);
            Populations = cell(1,K);
            GV          = cell(1,K);
            FrontNo     = cell(1,K);
            CrowdDis    = cell(1,K);
            index       = randperm(floor(Problem.N/K)*K);
            temp        = reshape(index,K,floor(Problem.N/K));
            for i = 1 : K
                Populations{i} = Population(temp(i,:));
                Masks{i}       = Mask(temp(i,:),:);
                Decs{i}        = Dec(temp(i,:),:);
                [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},length(Populations{i}));
                GV{i}          = UpdateGV(zeros(1,Problem.D),Masks{i},FrontNo{i});
            end
            endingFlag = 0;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                LastPopulation = [Populations{:}];              
                allSolutions   = [Populations{:}];
                ref            = max(allSolutions.objs);
                perHV          = zeros(1, K);
                for i = 1 : K
                    perHV(i) = HV(Populations{i}, ref);
                end
                minPopNum = max(2,floor(Problem.N/K*(Problem.FE/Problem.maxFE)));
                resources = allocateResourcesWithMinimum(perHV, Problem.N, minPopNum);
                for i = 1 : K                
                    GV{i}            = UpdateGV(GV{i},Masks{i},FrontNo{i});
                    Mating           = TournamentSelection(2,2*resources(i),FrontNo{i},-CrowdDis{i});
                    [OffDec,OffMask] = Operator(Problem,Decs{i}(Mating,:),Masks{i}(Mating,:),GV{i});
                    Offspring        = Problem.Evaluation(OffDec.*OffMask);
                    Populations{i}   = [Populations{i},Offspring];
                    Decs{i}          = [Decs{i};OffDec];
                    Masks{i}         = [Masks{i};OffMask];
                    if i > 1 && Problem.FE/Problem.maxFE > 1
                        for j = 1 : i-1
                            [~,fs(j)] = min(mean(Populations{j}.objs,2));
                        end
                        R = zeros(1,Problem.D);
                        for j = 1 : i-1
                            R = R + Masks{j}(fs(j),:);
                        end
                        R(R>0) = 1;
                        dis = sum(repmat(R,length(Populations{i}),1)&Masks{i},2);
                        [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},floor(Problem.N/K),dis);
                    else
                        [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},floor(Problem.N/K));
                    end
                end
                 
                CurrentPopulation = [Populations{:}];              
                globalSS          = simility(sum(LastPopulation.decs~=0)>mean(sum(LastPopulation.decs~=0)),sum(CurrentPopulation.decs~=0)>mean(sum(CurrentPopulation.decs~=0)));
                if globalSS >= 1
                    divisionFlag = 1;
                else
                    divisionFlag = 0;
                end
               
                if mod(ceil(Problem.FE/Problem.N),20) == 0
                    [ss,index] = SubPopSimility(Populations,Masks);
                    if ss > 0.5
                        K = K-1;
                        i = index(1);
                        j = index(2);
                        [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection([Populations{i},Populations{j}],[Decs{i};Decs{j}],[Masks{i};Masks{j}],floor(Problem.N/K));
                        Populations(j) = [];
                        Decs(j)        = [];
                        Masks(j)       = [];
                        GV(j)          = [];
                        endingFlag     = endingFlag + 1;
                    elseif divisionFlag == 1
                        for i = 1 : K
                            [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},floor(Problem.N/(K+1)));
                        end
                        K    = K + 1;
                        Dec  = unifrnd(repmat(Problem.lower,floor(Problem.N/K),1),repmat(Problem.upper,floor(Problem.N/K),1));
                        Mask = zeros(floor(Problem.N/K),Problem.D);
                        F    = zeros(1,Problem.D);
                        for i = 1: K-1
                            F = F + GV{i};
                        end
                        for i = 1 : floor(Problem.N/K)
                            Mask(i,TournamentSelection(2,floor(rand*Problem.D),F)) = 1;
                        end
                        Populations{K} = Problem.Evaluation(Dec.*Mask);
                        Masks{K}       = Mask;
                        Decs{K}        = Dec;
                        GV{K}          = zeros(1,Problem.D);
                        [Populations{K},Decs{K},Masks{K},FrontNo{K},CrowdDis{K}] = EnvironmentalSelection(Populations{K},Decs{K},Masks{K},length(Populations{K}));
                        GV{K}          = UpdateGV(zeros(1,Problem.D),Masks{K},FrontNo{K});
                    end
                end
                Population = [Populations{:}];
            end
        end
    end
end


function resources = allocateResourcesWithMinimum(perHV, N, minResources)    
    M = length(perHV); % Number of subpopulations
    reserved    = minResources * ones(1, M);
    remaining_N = N - sum(reserved);
    % Check if the total resources are sufficient
    if remaining_N < 0
        error('Total resources are insufficient to guarantee the minimum allocation.');
    end
    epsilon = 1e-6; % Small value to avoid division by zero
    normHV  = (perHV - min(perHV)) / (max(perHV) - min(perHV) + epsilon);   
    invHV   = 1 - normHV;   
    weights = invHV / sum(invHV);
    remainingResources = round(remaining_N * weights);
    diff = remaining_N - sum(remainingResources);
    if diff ~= 0
        [~, idx] = max(weights); % Assign the difference to the largest weight
        remainingResources(idx) = remainingResources(idx) + diff;
    end
    resources = reserved + remainingResources;
end

function s = simility(subPop1,subPop2)
	s = sum(subPop1&subPop2)/min(sum(subPop1),sum(subPop2));
end