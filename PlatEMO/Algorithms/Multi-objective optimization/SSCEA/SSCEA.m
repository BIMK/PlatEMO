classdef SSCEA < ALGORITHM
% <multi/many> <real/integer>
% Subspace segmentation based co-evolutionary algorithm
% nSel ---  5 --- Number of selected solutions for decision variable clustering
% nPer --- 50 --- Number of perturbations on each solution for decision variable clustering

%------------------------------- Reference --------------------------------
% G. Liu, Z. Pei, N. Liu, and Y. Tian, Subspace segmentation based
% co-evolutionary algorithm for balancing convergence and diversity in
% many-objective optimization, Swarm and Evolutionary Computation, 2023.
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
            [nSel,nPer] = Algorithm.ParameterSet(5,50);
            lb = Problem.N/10;
            ub = Problem.N;
            beta = 1;
            CAsize = floor(lb+(ub-lb)*Problem.FE^2*beta/Problem.maxFE^2);

            %% Generate random population
            Population = Problem.Initialization();
            CA = UpdateCA([],Population,CAsize);
            DA = UpdateDA([],Population,Problem.N);
            
            %% Detect the group of each distance variable
            [DV,CV] = VariableClustering(Problem,Population,nSel,nPer);
            
            %% Optimization
            while Algorithm.NotTerminated(DA)
                [~,D] = size(DA.decs);
                [ParentC,ParentM] = MatingSelection(CA,DA,Problem.N);
                if Problem.FE/Problem.maxFE < 0.5 || (Problem.FE/Problem.maxFE > 0.5 && rand <0.5)
                    OffDec = [ParentC,ParentM];
                    OffDec = OffDec.decs;
                    NewDec = [OperatorGA(Problem,ParentC.decs,{1,15,0,0});OperatorGA(Problem,ParentM.decs,{0,0,D/length(CV)/2,15})];
                    OffDec(:,CV) = NewDec(:,CV);
                    Offspring = Problem.Evaluation(OffDec);
                else
                    OffDec = [ParentC,ParentM];
                    OffDec = OffDec.decs;
                    NewDec = [OperatorGA(Problem,ParentC.decs,{1,15,0,0});OperatorGA(Problem,ParentM.decs,{0,0,D/length(DV)/2,15})];
                    OffDec(:,DV) = NewDec(:,DV);
                    Offspring = Problem.Evaluation(OffDec);
                end
                CAsize = floor(lb+(ub-lb)*Problem.FE^2*beta/Problem.maxFE^2);
                CA     = UpdateCA(CA,Offspring,CAsize);
                DA     = UpdateDA(DA,Offspring,Problem.N);
            end
        end
    end
end