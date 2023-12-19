classdef MGCEA < ALGORITHM
% <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Multi-granularity clustering based evolutionary algorithm

%------------------------------- Reference --------------------------------
% Y. Tian, S. Shao, G. Xie, and Y. Jin, A multi-granularity clustering
% based evolutionary algorithm for large-scale sparse multi-objective
% optimization, Swarm and Evolutionary Computation, 2023.
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
            %% Population initialization      
            [~,Fitness,SparseRate,TDec,TMask,TempPop] = FitnessCal(Problem);
            [Population,Dec,Mask,FitnessSpea2] = EnvironmentalSelection(TempPop,TDec,TMask,Problem.N);            
            NearStage = ceil(Problem.FE/(Problem.maxFE/10));
            [FitnessLayer,LayerMax] = UpdateLayer(SparseRate,NearStage,Fitness,Problem,[]);

            %% Optimization           
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,2*Problem.N,FitnessSpea2);
                [NearStage,Fitness,FitnessLayer,LayerMax] = ControlStage(SparseRate,NearStage,Mask,Dec,Fitness,FitnessLayer,LayerMax,Problem);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),FitnessLayer,LayerMax);
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FitnessSpea2] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
            end                       
        end
    end
end