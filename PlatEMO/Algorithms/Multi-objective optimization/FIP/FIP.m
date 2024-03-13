classdef FIP < ALGORITHM
% <multi> <real/integer>
% Regularity model-based multiobjective estimation of distribution
% algorithm

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, and Y. Jin, RM-MEDA: A regularity model-based
% multiobjective estimation of distribution algorithm, IEEE Transactions on
% Evolutionary Computation, 2008, 12(1): 41-63.
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
            Population = Problem.Initialization();
            obj1=[];
            obj2=[];
            %% Optimization
            while Algorithm.NotTerminated(Population)
                PopDec = Population.decs;
                obj1 = Problem.CalObj(ones(size(PopDec(1,:))));
                if isequal(obj1,obj2) ||isempty(obj2)
                    disp("1")
                else % 目标函数发生变化
                    FIP_IPG(Problem,PopDec);
                end
                obj2 = Problem.CalObj(ones(size(PopDec(1,:))));

                Offspring  = Operator(Problem,Population);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end