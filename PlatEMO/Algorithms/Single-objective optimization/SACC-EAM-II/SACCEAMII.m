classdef SACCEAMII < ALGORITHM
% <single> <real/integer> <expensive>
% Surrogate-assisted cooperative co-evolutionary algorithm of Minamo
% s --- 50 --- Number of subcomponents of separable variables

%------------------------------- Reference --------------------------------
% J. Blanchard, C. Beauthier, T. Carletti, A surrogate-assisted cooperative
% co-evolutionary algorithm using recursive differential grouping as
% decomposition strategy, Proceedings of the IEEE Congress on Evolutionary
% Computation, 2019: 689-696.
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
            s = Algorithm.ParameterSet(50);
            
            %% Separable decision variable grouping
            BU = Problem.lower;
            BD = Problem.upper;
            D  = Problem.D;
            [Population,Groups,K] = spiltVariables(Problem,BU,BD,s);
            
            %% Initialize each component
            % Start point
            GlobalIndi  = rand(1,D).*(BU-BD) + BD;
            GlobalObj   = Problem.Evaluation(GlobalIndi);
            Population  = [Population,GlobalObj];
            Fmin        = GlobalObj.objs;
            GlobalCandi = GlobalIndi;
            % Otain best solution of each component
            database = cell(1,K);
            for i = 1 : K
                select = Groups == i;
                subD   = sum(select);
                subN   = subD + 1;
                RepGlobal = repmat(GlobalIndi,subN,1);
                Decs   = rand(subN,subD).*repmat(BU(select)-BD(select),subN,1) + repmat(BD(select),subN,1);
                RepGlobal(:,select) = Decs;
                news   = Problem.Evaluation(RepGlobal);
                database{i} = news;
                Population  = [Population,news];
                [~,best]    = min(news.objs);
                BestDec     = news(best).decs;
                GlobalCandi(select) = BestDec(select);
            end
            % Evaluate the global candidate
            NewCandi   = Problem.Evaluation(GlobalCandi);
            CandiObj   = NewCandi.objs;
            Population = [Population,NewCandi];
            if CandiObj < Fmin
                Fmin = CandiObj;
                GlobalIndi = GlobalCandi;
            end
            
            %% Optimize
            while Algorithm.NotTerminated(Population)
                for i = 1 : K
                    tGlobalIndi = GlobalIndi;
                    select = Groups == i;
                    Decs   = database{i}.decs;
                    Objs   = database{i}.objs;
                    [model,~] = rbf_build(Decs(:,select),Objs);
                    PopSize   = 10*sum(select);
                    ProposedPoint       = GAOptimize(model,PopSize,Decs(:,select),BU(select),BD(select));
                    tGlobalIndi(select) = ProposedPoint;
                    newPropose    = Problem.Evaluation(tGlobalIndi);
                    database{i}   = [database{i},newPropose];
                    if newPropose.objs < CandiObj
                        decs = newPropose.decs;
                        GlobalCandi(select) = decs(select);
                        CandiObj = newPropose.objs;
                    end
                    Population = [Population,newPropose];
                end
                new = Problem.Evaluation(GlobalCandi);
                Population = [Population,new];
                if new.objs < Fmin
                    Fmin = new.objs;
                    GlobalIndi = GlobalCandi;
                end
            end
        end
    end
end