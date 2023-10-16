classdef NBLEA < ALGORITHM
% <multi> <real> <constrained/none> <bilevel>
% Nested bilevel evolutionary algorithm

%------------------------------- Reference --------------------------------
% A. Sinha, P. Malo, K. Deb, Test problem construction for single-objective 
% bilevel optimization, Evolutionary Computation, 2014, 22(3): 439-477.
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
            for i = 1 : size(ulPopDec,1)
               llPopDec(i,:) = llSearch(Problem,ulPopDec(i,:),[]);
            end
            Population = Problem.Evaluation([ulPopDec,llPopDec]); 
            
            %% UL Optimization
            while Algorithm.NotTerminated(Population)
                % Upper level optimization
                MatingPool = TournamentSelection(2,3,CalFitness(Problem.C,Population));
                ParentDec  = Population(MatingPool).decs;
                ulOffDec   = OperatorPCX(ParentDec(:,1:Problem.DU),Problem.lower(1:Problem.DU),Problem.upper(1:Problem.DU));
                % Lower level optimization
                [~,closest] = min(pdist2(ulOffDec,ParentDec(:,1:Problem.DU)),[],2);
                for i = 1 : size(ulOffDec,1)
                    llOffDec(i,:) = llSearch(Problem,ulOffDec(i,:),ParentDec(closest(i),Problem.DU+1:end));  
                end
                Offspring = Problem.Evaluation([ulOffDec,llOffDec]); 
                % Environment selection for Upper Population
                Population = EnvironmentalSelection(Problem,Population,Offspring);
            end
        end
	end
end