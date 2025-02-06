classdef NSBiDiCo < ALGORITHM
% <2023> <multi> <real/integer/label/binary/permutation> <constrained>
% Non-dominated sorting bidirectional differential coevolution algorithm

%------------------------------- Reference --------------------------------
% C. S. R. Mendes, A. F. R. Araujo, and L. R. C. Farias. Non-dominated
% sorting bidirectional differential coevolution. Proceedings of the IEEE
% International Conference on Systems, Mans and Cybernetics, 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias (email: lrcf@cin.ufpe.br)

	methods
		function main(Algorithm,Problem)
			%% Parameter setting
			[Cr,F,proM,disM] = Algorithm.ParameterSet(1,0.5,1,20);
			
            %% Generate random population
			Population = Problem.Initialization();
			ArcPop     = [];
			
            %% Optimization
			while Algorithm.NotTerminated(Population)		
				AllPop     = [Population,ArcPop];
				MatingPool = MatingSelection(Population,ArcPop,Problem.N);   
				Offspring  = OperatorDE(Problem, MatingPool(1:end), MatingPool(randi(Problem.N,1,Problem.N)), MatingPool(randi(Problem.N,1,Problem.N)), {Cr, F, proM, disM});
				ArcPop     = UpdateArc([AllPop, Offspring], Problem.N);				
				Population = EnvironmentalSelection([Population, Offspring], Problem.N);	
			end
        end
	end
end