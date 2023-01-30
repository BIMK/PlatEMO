classdef MOEARE < ALGORITHM  
% <multi> <real/integer/label/binary/permutation> <robust>
% Multi-objective evolutionary algorithm with robustness enhancement
% alpha --- 1.5 --- Parameter for deleting solutions from archive

%------------------------------- Reference --------------------------------
% Z. He, G. G. Yen, and J. Lv, Evolutionary multiobjective optimization
% with robustness enhancement, IEEE Transactions on Evolutionary
% Computation, 2020, 24(3): 494-507.
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
            alpha = Algorithm.ParameterSet(1.5);
            
            %% Generate random population
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N,Population.objs);
            z        = min(Population.objs,[],1);
            Archive  = Population;
            [~,ArcW] = min(pdist2(Archive.objs-z,W,'cosine'),[],2);
            ArcW     = num2cell(ArcW);
            ArcSP    = num2cell(sum(Archive.objs,2));

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Optimization
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                z          = min([z;Offspring.objs],[],1);
                [PopObjX,OffObjX,ArcObjX]             = Distrubance(Problem,Population.decs,Offspring.decs,Archive.decs);
                [Population,FrontNo,CrowdDis,PopObjX] = EnvironmentalSelection([Population,Offspring],Problem.N,[PopObjX;OffObjX]);
                % Update archive
                [SOI,SOIObjX] = SelectSOI(Population,PopObjX,20,z);
                [Archive,ArcW,ArcSP] = UpdateArchive(Archive,ArcObjX,SOI,SOIObjX,alpha,z,W,ArcW,ArcSP);
                % Final robust solution selection
                if Problem.FE >= Problem.maxFE
                    Population = FinalSelection(Archive,W,ArcW,ArcSP);
                end
            end
        end
    end
end