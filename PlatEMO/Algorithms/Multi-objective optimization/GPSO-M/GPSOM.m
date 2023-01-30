classdef GPSOM < ALGORITHM
% <multi/many> <real> <large/none> <constrained/none>
% Gradient based particle swarm optimization algorithm (for multi-objective optimization)
% popsize --- 20 --- Population size of single run

%------------------------------- Reference --------------------------------
% M. M. Noel, A new gradient based particle swarm optimization algorithm
% for accurate computation of global minimum, Applied Soft Computing, 2012,
% 12: 353-359.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            popsize = Algorithm.ParameterSet(20);

            %% Generate random population 
            SubPops = cell(1,Problem.N);
            Pbests  = cell(1,Problem.N);
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            SubPops{1}    = Problem.Initialization(popsize);
            Pbests(1)     = SubPops(1);
            FrontNo       = NDSort(SubPops{1}.objs,SubPops{1}.cons,1);
            index = FrontNo == 1;
            FrontNoPops = SubPops{1}(index);
            Gbest(1)    = LocalSearch(Problem,FrontNoPops(1),W(1,:));
            for i = 2 : Problem.N
                SubPops{i}  =  Problem.Initialization(popsize);
                Pbests(i)   = SubPops(i);
                [FrontNo,~] = NDSort(SubPops{i}.objs,SubPops{i}.cons,1);
                index = FrontNo == 1;
                FrontNoPops = SubPops{1}(index);
                Gbest(i)    = LocalSearch(Problem,FrontNoPops(1),W(i,:));
            end
            
            %% Optimization
            while Algorithm.NotTerminated(Gbest)
                for i = 1 : Problem.N
                    for j = 1 : popsize
                        SubPops{i}(j) = OperatorPSO(Problem,SubPops{i}(j),Pbests{i}(j),Gbest(i));
                        Pbests{i}(j)  = UpdatePbest(Pbests{i}(j),SubPops{i}(j));
                    end
                    FrontNo     = NDSort(SubPops{i}.objs,SubPops{i}.cons,1);
                    index       = FrontNo==1;
                    FrontNoPops = SubPops{i}(index);
                    Gbest(i)    = LocalSearch(Problem,FrontNoPops(1),W(i,:));
                end
            end
        end
	end
end