classdef RPEA < ALGORITHM
% <many> <real/integer/label/binary/permutation>
% Reference points-based evolutionary algorithm
% alpha --- 0.4 --- Ratio of individuals being used to generate reference points
% delta --- 0.1 --- Parameter determining the difference between the reference point and the individuals

%------------------------------- Reference --------------------------------
% Y. Liu, D. Gong, X. Sun, and Y. Zhang, Many-objective evolutionary
% optimization based on reference points, Applied Soft Computing, 2017,
% 50: 344-355.
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
            [alpha,delta] = Algorithm.ParameterSet(0.4,0.1);

            %% Generate random population
            Population = Problem.Initialization();
            R          = GenerateRefPoints(Population,delta*(max(Population.objs,[],1)-min(Population.objs,[],1)),alpha,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,min(TchebychevDistance(Population.objs,R),[],2));
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                R          = GenerateRefPoints([Population,Offspring],delta*(max(Population.objs,[],1)-min(Population.objs,[],1)),alpha,Problem.N);
                Population = EnvironmentalSelection([Population,Offspring],R,Problem.N);
            end
        end
    end
end