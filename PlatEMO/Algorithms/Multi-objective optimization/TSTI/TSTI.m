classdef TSTI < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Two-stage evolutionary algorithm with three indicators
% epsilon --- 0.05 --- parameter for calculating frank
% row     --- 1.01 --- parameter for adjusting epsilon

%------------------------------- Reference --------------------------------
% J. Dong, W. Gong, F. Ming, and L. Wang, A two-stage evolutionary
% algorithm based on three indicators for constrained multi-objective
% optimization, Expert Systems with Applications, 2022, 195: 116499.
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
			[Epsilon0,row] = Algorithm.ParameterSet(0.05,1.01);
			Population     = Problem.Initialization();
			[fpr,fcd]      = Estimation(Population.objs,1/Problem.N^(1/Problem.M));
			fcv            = Calculate_fcv(Population); 
			Epsilon        = Epsilon0;
			PopObj_1       = [fpr,fcd]; 
			[fm,~]         = NDSort(PopObj_1,Problem.N);
			PopObj         = [fm' + Epsilon * fcv,fcv];
			[frank,~]      = NDSort(PopObj,Problem.N);
			fitness        = frank' + fcv./(fcv+1);

			%% Optimization
			while Algorithm.NotTerminated(Population)
				if Problem.FE<=0.4*Problem.maxFE
					MatingPool = TournamentSelection(2,Problem.N,fitness);
					Offspring  = OperatorGA(Problem,Population(MatingPool));
					[fpr,fcd]  = Estimation(Offspring.objs,1/Problem.N^(1/Problem.M));
					fcv        = Calculate_fcv(Offspring); 
					OffObj_1   = [fpr,fcd]; 
					[fm,~]     = NDSort(OffObj_1,Problem.N);
					OffObj     = [fm' + Epsilon * fcv,fcv];
					[Population,fitness] = EnvironmentalSelectionStageI([Population,Offspring],PopObj,OffObj,Problem.N);
					[fpr,fcd]  = Estimation(Population.objs,1/Problem.N^(1/Problem.M));
					fcv        = Calculate_fcv(Population);
					PopObj_1   = [fpr,fcd]; 
					[fm,~]     = NDSort(PopObj_1,Problem.N);
					PopObj     = [fm' + Epsilon * fcv,fcv];
					Epsilon    = row * Epsilon;
				else
					[Population,fit2] = EnvironmentalSelectionStageII(Population,Problem.N);
					MatingPool        = TournamentSelection(2,Problem.N,fit2);
					Offspring         = OperatorGA(Problem,Population(MatingPool));
					[Population,fit2] = EnvironmentalSelectionStageII([Population,Offspring],Problem.N);
                end	
			end
		end
	end
end