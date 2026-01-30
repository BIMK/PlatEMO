classdef AMGPSL < ALGORITHM
% <2025> <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% Adaptive multi-granular Pareto-optimal subspace learning

%------------------------------- Reference --------------------------------
% C. Sun, Y. Tian, S. Shao, S. Yang, and X. Zhang. An adaptive multi-
% granular Pareto-optimal subspace learning algorithm for sparse large-
% scale multi-objective optimization. Proceedings of the IEEE Congress on
% Evolutionary Computation, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
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
            [Population,Dec,Mask,FitnessSpea2]        = EnvironmentalSelection(TempPop,TDec,TMask,Problem.N);
            
            %% Initialize stage parameters
            NearStage               = ceil(Problem.FE/(Problem.maxFE/10));
            [FitnessLayer,LayerMax] = UpdateLayer(SparseRate,NearStage,Fitness,Problem,[]);
            rho = 0.5;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Model training phase
                Site = rho > rand(1,ceil(Problem.N/2));
                if any(Site)
                    [rbm,dae,allZero,allOne] = ModelTraining(Mask,Dec,any(Problem.encoding~=4),Problem);
                else
                    [rbm,dae,allZero,allOne] = deal([]);
                end
                
                % Mating selection & Update layers
                MatingPool = TournamentSelection(2,2*Problem.N,FitnessSpea2);
                [NearStage,Fitness,FitnessLayer,LayerMax] = ControlStage(SparseRate,NearStage,Mask,Dec,Fitness,FitnessLayer,LayerMax,Problem);
                
                % Generate offspring
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),rbm,dae,Site,allZero,allOne,FitnessLayer,LayerMax);
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                
                % Environmental selection
                [Population,Dec,Mask,FitnessSpea2] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
                
                % Update rho adaptively
                success = false(1,length([Population,Offspring]));
                success(1:length(Population)) = true;
                s1  = sum(success(length(Population)+1:end))/2;
                s2  = sum(~success(length(Population)+1:end))/2;
                rho = (rho + min(max((s1+1e-6)/(s1+s2+1e-6),0.1),0.9))/2;
            end
        end
    end
end