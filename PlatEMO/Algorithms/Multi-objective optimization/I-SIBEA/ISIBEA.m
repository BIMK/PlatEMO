classdef ISIBEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% Interactive simple indicator-based evolutionary algorithm
% Point --- --- Preferred point

%------------------------------- Reference --------------------------------
% T. Chugh, K. Sindhya, J. Hakanen, and K. Miettinen, An interactive simple
% indicator-based evolutionary algorithm (I-SIBEA) for multiobjective
% optimization problems, Proceedings of the International Conference on
% Evolutionary Multi-Criterion Optimization, 2015, 277-291.
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
            Point = Algorithm.ParameterSet(ones(1,Problem.M));

            %% Generate random population
            Population = Problem.Initialization();
            FrontNo    = NDSort(Evaluate(Population.objs,Point),inf);
            WHVLoss    = CalWHVLoss(Population.objs,FrontNo);
            wz = [];    % Weight distribution function value
            AA = [];    % Preferred set
            RA = [];    % Non-preferred set

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-WHVLoss);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,FrontNo,WHVLoss] = EnvironmentalSelection([Population,Offspring],Problem.N,wz,AA,RA);
                if ~mod(ceil(Problem.FE/Problem.N),ceil(ceil(Problem.maxFE/Problem.N)/4))
                    [wz,AA,RA] = Interaction(Population.objs,Point);
                end
            end
        end
     end
end

function PopObj = Evaluate(PopObj,Point)
% g-dominance based function evaluation

    Point = repmat(Point,size(PopObj,1),1);
    Flag  = all(PopObj<=Point,2) | all(PopObj>=Point,2);
    Flag  = repmat(Flag,1,size(PopObj,2));
    PopObj(~Flag) = PopObj(~Flag) + 1e10;
end