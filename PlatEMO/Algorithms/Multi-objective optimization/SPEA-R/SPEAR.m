classdef SPEAR < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Strength Pareto evolutionary algorithm based on reference direction

%------------------------------- Reference --------------------------------
% S. Jiang and S. Yang, A strength Pareto evolutionary algorithm based on
% reference direction for multiobjective and many-objective optimization,
% IEEE Transactions on Evolutionary Computation, 2017, 21(3): 329-346.
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
            %% Generate the reference directions (general approach)
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            % Largest acute angle between two neighbouring reference directions
            cosine = 1 - pdist2(W,W,'cosine');
            cosine(logical(eye(length(cosine)))) = 0;
            theta  = max(min(acos(cosine),[],2));

            %% Generate random population
            Population = Problem.Initialization();

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs,20);
                Offspring  = OperatorGAhalf(Problem,Population([1:Problem.N,MatingPool]));
                QObj       = ObjectiveNormalization([Population,Offspring]);
                [Ei,Angle] = Associate(QObj,W);
                FV         = FitnessAssignment(Ei,QObj,Angle,theta);
                Population = EnvironmentalSelection(Problem,[Population,Offspring],Ei,FV);
            end
        end
    end
end