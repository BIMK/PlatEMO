classdef CCGDE3 < ALGORITHM
% <multi> <real/integer> <large/none>
% Cooperative coevolution generalized differential evolution 3

%------------------------------- Reference --------------------------------
% L. M. Antonio and C. A. Coello Coello, Use of cooperative coevolution for
% solving large scale multiobjective optimization problems, Proceedings of
% the IEEE Congress on Evolutionary Computation, 2013, 2758-2765.
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
            %% Random grouping
            NumEsp = 2; % the number of subPopulations
            varsPerGroup = floor(Problem.D/NumEsp);
            Index = [];
            numSPop = ceil(0.8*Problem.N);  % number of individuals in each subPopulation
            for i = 1:NumEsp-1
                Index = [Index,ones(1,varsPerGroup).*i];
            end
            Index = [Index, ones(1,Problem.D-size(Index,2)).*NumEsp];
            Index = Index(randperm(length(Index))); % randomly grouping
            
            %% Generate random population
            Gmax        = 1;	% the times of each subpopulation
            Dec         = unifrnd(repmat(Problem.lower,numSPop,1),repmat(Problem.upper,numSPop,1));
            subDec1     = Dec(:,Index==1);
            subDec2     = Dec(:,Index==2);
            Population1 = GetInd(Problem,subDec1,subDec2,Index,numSPop,1);
            Population2 = GetInd(Problem,subDec2,subDec1,Index,numSPop,2);
            Population  = EnvironmentalSelection([Population1,Population2],Problem.N);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                FrontNo1  = NDSort(Population1.objs,1);
                FrontNo2  = NDSort(Population2.objs,1);
                NDsubDec1 = subDec1(FrontNo1==1,:);
                NDsubDec2 = subDec2(FrontNo2==1,:);
                for j = 1 : NumEsp
                    for k = 1 : Gmax
                        if j == 1
                        	OffDec1     = CCDE(subDec1,subDec1(randi(numSPop,1,numSPop),:),subDec1(randi(numSPop,1,numSPop),:),Problem.lower(Index==j),Problem.upper(Index==j));
                        	Offspring   = GetInd(Problem,OffDec1,NDsubDec2,Index,numSPop,j);
                        	Population1 = GDE3_EnvironmentalSelection(Population1,Offspring,numSPop);
                        	Dec1        = Population1.decs;
                        	subDec1     = Dec1(:,Index==1);
                        else
                        	OffDec2     = CCDE(subDec2,subDec2(randi(numSPop,1,numSPop),:),subDec2(randi(numSPop,1,numSPop),:),Problem.lower(Index==j),Problem.upper(Index==j));
                        	Offspring   = GetInd(Problem,OffDec2,NDsubDec1,Index,numSPop,j);
                        	Population2 = GDE3_EnvironmentalSelection(Population2,Offspring,numSPop);
                        	Dec2        = Population1.decs;
                        	subDec2     = Dec2(:,Index==2);
                        end
                        Population = EnvironmentalSelection([Population1,Population2],Problem.N);
                    end
                end
            end
        end
    end
end

function Population = GetInd(Problem,subDec1,subDec2,Index,numSPop,j)
    Dec = zeros(numSPop,Problem.D);
    if j == 1
       Dec(:,Index==1) = subDec1;
       Dec(:,Index==2) = subDec2(randi(end,numSPop,1),:);
    else
       Dec(:,Index==2) = subDec1;
       Dec(:,Index==1) = subDec2(randi(end,numSPop,1),:);           
    end
    Population = Problem.Evaluation(Dec);
end