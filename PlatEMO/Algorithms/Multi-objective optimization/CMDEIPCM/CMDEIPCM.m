classdef CMDEIPCM < ALGORITHM
% <2022> <multi> <real/integer> <large/none> <constrained>
% Constrained multiobjective differential evolution algorithm with an infeasible proportion control mechanism

%------------------------------- Reference --------------------------------
% J. Liang, X. Ban, K. Yu, K. Qiao, and B. Qu. Constrained multiobjective
% differential evolution algorithm with infeasible-proportion control
% mechanism. Knowledge-Based Systems, 2022, 250: 109105.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();    
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Fitness2    = CalFitness(Population2.objs);
    
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                Offspring1 = DEgenerator(Population1,Population2,Problem,Fitness1);
                Offspring2 = DEgenerator(Population2,Population1,Problem,Fitness2);
                [Population1,Fitness1] = EnvironmentalSelection1([Population1,Offspring1,Offspring2],Problem.N,true);
                pinfea     = 0.5*(1-cos((1-Problem.FE./Problem.maxFE)*pi));
                Population = [Population2,Offspring1,Offspring2];
                Obj        = Population.objs;
                MaxObj     = max(Obj);
                cons  = Population.cons;
                cons(cons<=0) = 0;
                CV    = sum(cons,2);
                Infea = find(CV > 0);
                p_in  = size(Infea,1)/size(Population.decs,1);
                if p_in > pinfea
                    Infea_Sec = randperm(size(Infea,1),floor(size(Population.decs,1)*(p_in - pinfea)));
                    Obj(Infea(Infea_Sec),:) = Obj(Infea(Infea_Sec),:) + MaxObj;
                end
                [Population2,Fitness2] = EnvironmentalSelection2(Population,Obj,Problem.N,false);
            end
        end
    end
end