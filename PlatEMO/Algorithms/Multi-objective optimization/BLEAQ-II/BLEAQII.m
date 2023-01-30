classdef BLEAQII < ALGORITHM
% <multi> <real> <constrained/none> <bilevel>
% Bilevel evolutionary algorithm based on quadratic approximations II

%------------------------------- Reference --------------------------------
% A. Sinha, Z. Lu, K. Deb, and P. Malo, Bilevel optimization based on
% iterative approximation of mappings, Journal of Heuristics, 2020, 26:
% 151-185.
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
            %% Generate random population
            ulPopDec = unifrnd(repmat(Problem.lower(1:Problem.DU),Problem.N,1),repmat(Problem.upper(1:Problem.DU),Problem.N,1));
            % Doing lower level search for population member
            for i = 1 : size(ulPopDec,1)
               [llPopDec(i,:),tag.ulPop(i)] = llSearch(Problem,ulPopDec(i,:),[],[],[]);
            end
            llLocalSearchPopSize = max([4 Problem.N/10]);
            [~,centroids] = kmeans(llPopDec,llLocalSearchPopSize,'EmptyAction','singleton');
            for i = 1 : size(ulPopDec,1)
               [llPopDec(i,:),tag.ulPop(i)] = llSearch(Problem,ulPopDec(i,:),[],[], centroids);
            end
            Population           = Problem.Evaluation([ulPopDec,llPopDec]);
            archive              = struct('tag1',[],'tag0',[]);
            archive.tag1 = archiveUpdate(Problem, archive.tag1, Population);
            alphaStoppingInitial = sum(var([ulPopDec(tag.ulPop==1,:), llPopDec(tag.ulPop==1,:)]))/Problem.D;
            maxError             = 1e-4;
            StoppingCriteria     = 0;
            gen = 0;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen = gen+1;
                % Upper level optimization
                MatingPool = TournamentSelection(2,3,CalFitness(Problem.C,Population));
                ParentDec  = Population(MatingPool).decs;
                ulOffDec   = OperatorPCX(ParentDec(:,1:Problem.DU),Problem.lower(1:Problem.DU),Problem.upper(1:Problem.DU));
                
                % Lower level optimization
                for i = 1 : size(ulOffDec,1)
                    if (sum(tag.ulPop==1) < Problem.N/2) || (length(archive.tag1) < (Problem.DU+1)*(Problem.DU+2)/2+3*(Problem.DU))
                        % Optimization
                        [~,closest]   = min(pdist2(ulOffDec,ParentDec(:,1:Problem.DU)),[],2);
                        [llOffDec(i,:), tag.Offsprings(i)] = llSearch(Problem,ulOffDec(i,:),ParentDec(closest(i),Problem.DU+1:end),[],[]);  
                    else
                        % Approximation
                        [psiMapping,phiMapping,lies]    = getMappings(Problem,ulOffDec(i,:),archive.tag1);
                        [llOffDec(i,:),sumMSE,validMSE] = getLowerLevelVariableFromMapping(ulOffDec(i,:),psiMapping,phiMapping,Problem,archive);
                        if lies==1 && validMSE<maxError
                            tag.Offsprings(i) = 1;
                        else
                            tag.Offsprings(i) = 0;
                        end
                    end
                end
                Offspring = Problem.Evaluation([ulOffDec,llOffDec]);
                
                % Environment selection for upper population
                Population       = EnvironmentalSelection(Problem,Population,Offspring);
                
                PopDec           = Population.decs;
                llmemberVariance = var(PopDec(:,Problem.DU+1:end));
                if sum(tag.Offsprings == 1) > 0 
                    archive.tag1 = archiveUpdate(Problem, archive.tag1, Offspring(tag.Offsprings==1));
                end
                if sum(tag.Offsprings == 0) > 0
                    archive.tag0 = archiveUpdate(Problem, archive.tag0, Offspring(tag.Offsprings==0));
                end
                
                % Local search in EliteIndiv
                Fitness    = CalFitness(Problem.C,Population);
                if sum(tag.ulPop==1)==0
                    [~,Index] = min(Fitness);
                else
                    [~,index]  = min(Fitness(tag.ulPop==1));
                    I          = find(tag.ulPop==1);
                    Index      = I(index);
                end
                initialIndivLS = Population(Index);
                PopDec        = Population.decs;
                llPopDec      = PopDec(:,Problem.DU+1:end); 
                alphaStopping = sum(var([ulPopDec(tag.ulPop==1,:),llPopDec(tag.ulPop==1,:)]))/Problem.D;
                alphaStopping = alphaStopping/alphaStoppingInitial;
                if alphaStopping < 1e-4
                    StoppingCriteria = 1;
                end
                [eliteIndivLS,llEliteIndivLS] = DoLocalSearch(Problem,initialIndivLS,archive,StoppingCriteria,gen);
                [llEliteIndivLS,           ~] = llSearch(Problem,eliteIndivLS,llEliteIndivLS,llmemberVariance,[]);
                eliteLS = Problem.Evaluation([eliteIndivLS,llEliteIndivLS]);
                % EnvironmentalSelection
                if CalFitness(Problem.C,eliteLS) > CalFitness(Problem.C,initialIndivLS)
                    eliteLS = initialIndivLS;
                end
                Population = EnvironmentalSelection(Problem,Population,eliteLS);
            end
        end
    end
end

function archive = archiveUpdate(Problem, archive, newData)

    archiveSize = ((Problem.DU+1)*(Problem.DU+2)/2+3*(Problem.DU))*10;
    archive = [archive, newData];
    
    if length(archive) > archiveSize
        archive(1) = [];
    end 
end