classdef PREA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Promising-region based EMO algorithm

%------------------------------- Reference --------------------------------
% J. Yuan, H. Liu, F. Gu, Q. Zhang, and Z. He, Investigating the properties
% of indicators and an evolutionary many-objective algorithm based on a
% promising region, IEEE Transactions on Evolutionary Computation, 2021,
% 25(1): 75-86.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiawei Yuan

    methods
        function main(Algorithm,Problem)
            %% Generate the random population
            Population = Problem.Initialization();
            Zmin       = min(Population.objs,[],1);

            %% shift the objective space to R+
            PopObj = Population.objs;
            PopObj = PopObj - repmat(Zmin,Problem.N,1) + 1e-6;

            %% calculate the ratio based indicator matrix
            IMatrix = ones(Problem.N,Problem.N); 
            for i=1:1:Problem.N
                Fi             = PopObj(i,:);
                % calculate ratio based indicator value of each individual
                Ir             = PopObj./repmat(Fi,Problem.N,1) - 1;  
                InvertIr       = repmat(Fi,Problem.N,1)./PopObj - 1;
                MaxIr          = max(Ir,[],2);
                MinIr          = max(InvertIr,[],2);
                DomInds        = find(MaxIr<=0);
                MaxIr(DomInds) = -MinIr(DomInds);
                IMatrix(i,:)   = MaxIr';
                IMatrix(i,i)   = Inf;
            end

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingStrategy(IMatrix);
                Offspring  = OperatorGAhalf(Problem,Population(MatingPool));
                Zmin       = min([Zmin;Offspring.objs],[],1);
                [Population,IMatrix] = PREA_Update([Population,Offspring],Problem.N,Zmin);
            end
        end
    end
end