classdef NRVMOEA < ALGORITHM
% <2024> <multi/many> <real/integer/label/binary/permutation>
% Adaptive normal reference vector-based multi- and many-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% Y. Hua, Q. Liu, and K. Hao. Adaptive normal vector guided evolutionary
% multi- and many-objective optimization. Complex & Intelligent Systems,
% 2024, 10: 3709-3726.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yicun Hua

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            EP         = [];

            Znadir      = max(Population.objs,[],1);
            Zmin        = min(Population.objs,[],1);
            Zmax        = max(Population.objs,[],1);
            scale       = Zmax-Zmin;
            Archive     = UpdateArchive(Population,[],Problem.N);
            extremPoint = ones(Problem.M,Problem.M).*10e30;

            PopObjn       = (Population.objs - Zmin )./scale;
            ArcObjn       = (Archive.objs - Zmin )./scale;
            Hyperplane_bp = [];

            [~,Hyperplane_bp] = vertmap(ArcObjn,PopObjn,Hyperplane_bp,Problem);
            Hyperplane        = Hyperplane_bp;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                Archive    = UpdateArchive(Population,Archive,Problem.N);
                Population = [Population Archive];
                MatingPool = randi(length(Population),1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                UniPop     = [Population,Offspring];
                PopObj     = UniPop.objs;
                [FrontNo,MaxFNo] = NDSort(PopObj,Problem.N);
                ParetoP    = find(FrontNo==MaxFNo);
                chosenPop  = UniPop(FrontNo<MaxFNo);
                id = find(FrontNo<MaxFNo);
                id = id(randperm(size(id,2),size(id,2)));
                otherPopObj = UniPop(id).objs;
                if size(ParetoP,2) < Problem.M
                    PopObj  = UniPop([ParetoP id(1:10-size(ParetoP,2))]).objs;
                    ParetoP = [ParetoP id(1:10-size(ParetoP,2))];
                else
                    PopObj = UniPop(ParetoP).objs;
                end
                Zmin = min(UniPop.objs,[],1);
                Zmax = max(PopObj,[],1);
                Zmin = min([Zmin;PopObj],[],1);
                [Znadir,extremPoint] = updateNadirPoint(Archive.objs,Zmin,Znadir,extremPoint);
                if ~mod(Problem.FE,ceil(0.1*Problem.maxFE))
                    scale = Zmax-Zmin;
                    scale(scale==0) = 10^(-5);
                    Hyperplane_bp   = Hyperplane;
                end
                scale(scale==0) = 10^(-5);
                PopObjn = (PopObj-repmat(Zmin,size(PopObj,1),1))./repmat(scale,size(PopObj,1),1);
                ArcObjn = (Archive.objs-repmat(Zmin,size(Archive.objs,1),1))./repmat(scale,size(Archive.objs,1),1);
                [mapPop,Hyperplane] = vertmap(ArcObjn,PopObjn,Hyperplane_bp,Problem);
                try
                    T = clusterdata(mapPop,'maxclust',Problem.N-length(chosenPop),'distance','euclidean','linkage','ward');
                catch e
                    T = clusterdata(PopObjn,'maxclust',Problem.N-length(chosenPop),'distance','euclidean','linkage','ward');
                end
                for c = 1 : Problem.N-length(chosenPop)
                    current = find(T == c);
                    pn      = length(current);
                    Ref     = sum(mapPop(current,:),1)/pn;
                    if pn > 1
                        d12 = zeros(pn,1);
                        for pc = 1 : pn
                            d1 = norm(Ref-mapPop(current(pc),:));
                            d2 = -(PopObjn(current(pc),:)*Hyperplane-1)./sqrt(sum(Hyperplane.^2));
                            d12(pc) = d1-d2;
                        end
                        [~,ct] = min(d12);
                        choose = current(ct);
                    else
                        choose = current;
                    end
                    EP = [EP,UniPop(ParetoP(choose)),chosenPop];
                    EP = unique(EP);
                end
                if length(EP) > Problem.N || Problem.FE >= 0.9*Problem.maxFE
                    [EPFNo,EPMaxFNo] = NDSort(EP.objs,Problem.N);
                    EP    = EP(EPFNo<=EPMaxFNo);
                    EPObj = EP.objs;
                    [~,Rank]   = sort(EPObj);
                    Extreme    = zeros(1,Problem.M);
                    Extreme(1) = Rank(1,1);
                    for j = 2 : length(Extreme)
                        k = 1;
                        Extreme(j) = Rank(k,j);
                        while ismember(Extreme(j),Extreme(1:j-1))
                            k = k+1;
                            Extreme(j) = Rank(k,j);
                        end
                    end
                    EPtemp      = EP(Extreme);
                    EP(Extreme) = [];
                    LEP         = length(EP);
                    Extreme     = [];
                    for d = 1 : (LEP-Problem.N+length(Extreme))
                        EPObj   = EP.objs;
                        dis     = pdist2((EPObj-Zmin)./scale,(EPObj-Zmin)./scale);
                        dis(logical(eye(length(dis)))) = 10^10;
                        mindis  = min(dis,[],2);
                        [~,del] = min(mindis);
                        EP(del) = [];
                    end
                    EP = [EP,EPtemp];
                end
                Population = EP;
                EP         = [];
                PopObjn    = (Population.objs - Zmin )./scale;
            end
        end
    end
end