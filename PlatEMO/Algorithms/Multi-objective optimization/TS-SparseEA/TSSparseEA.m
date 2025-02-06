classdef TSSparseEA < ALGORITHM
% <2022> <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Two-stage SparseEA
% r_eva  --- 0.1 --- The ratio of evaluations for the group optimization
% nGroup ---  50 --- The group size for the group optimization 

%------------------------------- Reference --------------------------------
% J. Jiang, F. Han, J. Wang, Q. Ling, H. Han, and Y. Wang. A two-stage
% evolutionary algorithm for large-scale sparse multiobjective optimization
% problems. Swarm and Evolutionary Computation, 2022, 72: 101093.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter settings
            [r_eva,nGroup] = Algorithm.ParameterSet(0.1,50);
            
            %% Initialization
            if Problem.encoding(1) == 4
                REAL = 0;
            else
                REAL = 1;
            end
            if REAL
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            else
                Dec = ones(Problem.N,Problem.D);
            end
            Mask = binornd(ones(Problem.N,Problem.D),0.5);
            Population = Problem.Evaluation(Dec.*Mask);
            
            %% Optimization
            [Population,Dec,Mask]    = BinaryGroupOptimization(Problem,Population,Dec,Mask,r_eva,nGroup,REAL);
            [~,~,~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N);
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),REAL);
                if REAL
                    OffDec = Match(OffDec,OffMask,Problem);
                end
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
            end
        end
    end
end