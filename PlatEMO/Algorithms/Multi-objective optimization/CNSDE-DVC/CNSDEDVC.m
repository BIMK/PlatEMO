classdef CNSDEDVC < ALGORITHM
% <multi> <real/integer> <robust>
% Constrained nondominated sorting differential evolution based on decision variable classification
% SN     ---     4 --- Number of perturbed solutions
% PN     ---     6 --- Number of perturbations
% TN     ---    15 --- Number of repeated times of perturbation
% theta  --- 0.001 --- Threshold for DVC operation
% eta    --- 0.001 --- Desired level of robustness

%------------------------------- Reference --------------------------------
% W. Du, W. Zhong, Y. Tang, W. Du, and Y. Jin, High-dimensional robust
% multi-objective optimization for order scheduling: A decision variable
% classification approach, IEEE Transactions on Industrial Informatics,
% 15(1): 293-304.
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
            [SN,PN,TN,theta,eta] = Algorithm.ParameterSet(4,6,15,0.001,0.001);

            %% Generate random population
            Population = Problem.Initialization();
            [HR,LR]    = DVC(Problem,Population,SN,PN,TN,theta);
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Problem,Population,Problem.N,false,eta);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if ~isempty(HR)
                    for subgen = 1 : 10
                        OffDec       = Population(TournamentSelection(2,end,FrontNo,-CrowdDis)).decs;
                        NewDec       = OperatorDE(Problem,Population.decs,Population(randi(end,1,end)).decs,Population(randi(end,1,end)).decs,{0.9,0.5,1,20});
                        OffDec(:,HR) = NewDec(:,HR);
                        Offspring    = Problem.Evaluation(OffDec);
                        [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Problem,[Population,Offspring],Problem.N,true,eta);
                    end
                end
                if ~isempty(LR)
                    for subgen = 1 : 2
                        OffDec       = Population(TournamentSelection(2,end,FrontNo,-CrowdDis)).decs;
                        NewDec       = OperatorDE(Problem,Population.decs,Population(randi(end,1,end)).decs,Population(randi(end,1,end)).decs,{0.9,0.5,1,20});
                        OffDec(:,LR) = NewDec(:,LR);
                        Offspring    = Problem.Evaluation(OffDec);
                        [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Problem,[Population,Offspring],Problem.N,false,eta);
                    end
                end
            end
        end
    end
end