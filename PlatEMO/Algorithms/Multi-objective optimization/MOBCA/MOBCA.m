classdef MOBCA < ALGORITHM
% <2024> <multi> <real/integer>
% Multi-objective besiege and conquer algorithm
% BCB       --- 0.2 --- Set BCB
% nSoldiers ---   3 --- Number of soldiers for each armies
% div       ---  10 --- Division number of grids

%------------------------------- Reference --------------------------------
% J. Jiang, J. Wu, J. Luo, X. Yang, and Z. Huang. MOBCA: multi-objective
% besiege and conquer algorithm. Biomimetics, 2024, 9: 316.
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
            %% Generate random population
            [BCB,nSoldiers,div] = Algorithm.ParameterSet(0.2,3,10);
            nArmies             = fix(Problem.N/nSoldiers);

            %% Generate random population
            Population          = Problem.Initialization(nArmies);
            Archive             = Population;
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                REP        = REPSelection(Archive.objs,nArmies,div);
                Offspring  = OperatorBCAGrid(Problem,Population,BCB,Archive(REP),nSoldiers,nArmies);
                Archive    = UpdateArchive([Archive,Offspring],Problem.N,div);
                Population = BCAUpdatePop(Archive,Population,nArmies,div);
            end
        end
    end
end