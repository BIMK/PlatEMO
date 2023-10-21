classdef MSCEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Multi-stage constrained multi-objective evolutionary algorithm
% cp --- 5 --- Decrease trend of the dynamic constraint boundary

%------------------------------- Reference --------------------------------
% Y. Zhang, Y. Tian, H. Jiang, X. Zhang, and Y. Jin, Design and analysis of
% helper-problem-assisted evolutionary algorithm for constrained 
% multiobjective optimization, Information Sciences, 2023, 648: 119547.
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
            %% Parameter setting
            cp = Algorithm.ParameterSet(5);
            %% Generate the random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();            
            %% Calculate the initial dynamic constraint boundary
            [~, nCon]               = size(Population1.cons);
            [initialE1, ~]          = max(max(0,Population1.cons), [], 1);
            [initialE2, ~]          = max(max(0,Population2.cons), [], 1);
            initialE1(initialE1==0) = 1;            
            initialE2(initialE2==0) = 1; 
            epsn1                   = initialE1;
            epsn2                   = initialE2;    
            CV2                     = sum(max(0,Population2.cons),2);            
            Fitness2                = CalFitness([CalSDE(Population2.objs)',CV2]);
            MaxCV2                  = zeros(ceil(Problem.maxFE/Problem.N),1);
            MaxCV2(1)               = max(CV2);
            ASC2                    = 0;
            arch                    = archive([Population1,Population2],Problem.N);
            %% Optimization
            while Algorithm.NotTerminated(Population1)                
                PopCon1   = max(0,Population1.cons);
                PopCon2   = max(0,Population2.cons);
                if sum(sum(PopCon1<=epsn1,2)==nCon) == length(Population1)
                    epsn1 = ReduceBoundary(initialE1,ceil(Problem.FE/Problem.N),ceil(Problem.maxFE/Problem.N)-1,cp);
                end
                CV2       = sum(PopCon2,2);
                MaxCV2(ceil(Problem.FE/Problem.N)) = max(CV2);
                if MaxCV2(ceil(Problem.FE/Problem.N))-MaxCV2(ceil((Problem.FE-Problem.N)/Problem.N))>0
                    ASC2  = ASC2 + 1;
                    epsn2 = ReduceBoundary(initialE2,ceil(Problem.FE/Problem.N)-ASC2,ceil(Problem.maxFE/Problem.N)-1,cp);
                elseif sum(sum(PopCon2<=epsn2,2)==nCon) == length(Population2)
                    ASC2  = 0;
                    epsn2 = ReduceBoundary(initialE2,ceil(Problem.FE/Problem.N),ceil(Problem.maxFE/Problem.N)-1,cp);
                end                                                        
                MatingPool1 = TournamentSelection(2,Problem.N,sum(max(0,Population1.cons-epsn1),2));
                MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                Offspring2  = OperatorGAhalf(Problem,Population2(MatingPool2));
                Population1            = EnvironmentalSelection1([Population1,Offspring1,Offspring2],Problem.N,epsn1);
                [Population2,Fitness2] = EnvironmentalSelection2([Population2,Offspring1,Offspring2],Problem.N,epsn2);                
                % Output the non-dominated and feasible solutions.
                arch = [arch,Population1,Population2];
                 [~, Unduplicated] = unique(arch.objs,'rows');
                arch = arch(Unduplicated);
                arch = archive(arch,Problem.N);
                if Problem.FE >= Problem.maxFE
                    Population1 = arch;
                end                        
            end
        end
    end
end