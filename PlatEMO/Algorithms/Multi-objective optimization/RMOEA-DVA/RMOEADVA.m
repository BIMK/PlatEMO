classdef RMOEADVA< ALGORITHM
% <multi> <real/integer> <robust>
% Robust multi-objective evolutionary algorithm with decision variable assortment
% nDVA  ---  50 --- Number of solutions for decision variable assortment
% theta --- 0.3 --- Threshold for decision variable assortment

%------------------------------- Reference --------------------------------
% J. Liu, Y. Liu, Y. Jin, and F. Li, A decision variable assortment-based
% evolutionary algorithm for dominance robust multiobjective optimization,
% IEEE Transactions on Systems, Man, and Cybernetics: Systems, 2022, 52(5):
% 3360-3375.
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
            [nDVA,theta] = Algorithm.ParameterSet(50,0.3);
            
            %% Generate random population
            Population = Problem.Initialization();
            [HR,LR]    = DVA(Problem,Population,nDVA,theta);
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Problem,Population,Problem.N,false);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if ~isempty(LR)
                    OffDec       = Population(TournamentSelection(2,end,FrontNo,-CrowdDis)).decs;
                    NewDec       = OperatorGA(Problem,Population(randi(end,1,end)).decs);
                    OffDec(:,LR) = NewDec(:,LR);
                    Offspring    = Problem.Evaluation(OffDec);
                    [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Problem,[Population,Offspring],Problem.N,false);
                end
                if ~isempty(HR)
                    OffDec       = Population(TournamentSelection(2,end,FrontNo,-CrowdDis)).decs;
                    NewDec       = OperatorGA(Problem,Population(randi(end,1,end)).decs);
                    OffDec(:,HR) = NewDec(:,HR);
                    Offspring    = Problem.Evaluation(OffDec);
                    [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Problem,[Population,Offspring],Problem.N,true);
                end
            end
        end
    end
end