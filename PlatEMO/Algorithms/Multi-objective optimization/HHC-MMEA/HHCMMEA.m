classdef HHCMMEA < ALGORITHM
% <multi> <real> <large/none> <multimodal> <sparse>
% Hybrid hierarchical clustering based multi-modal multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% Z. Ding, L. Cao, L. Chen, D. Sun, X. Zhang, and Z. Tao, Large-scale
% multimodal multiobjective evolutionary optimization based on hybrid
% hierarchical clustering, Knowledge-Based Systems, 2023, 266: 110398.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%This function is written by Lve Cao
%Fitness vector is the global guide vector and GV is the local guide vector
%in the original paper;

    methods
        function main(Algorithm,Problem)
            %% Calculate the fitness of each decision variable
            Fitness = zeros(1,Problem.D);
            REAL    = ~strcmp(Problem.encoding,'binary');
            for i = 1 : 1+4*REAL
                if REAL
                    Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                else
                    Dec = ones(Problem.D,Problem.D);
                end
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                Fitness    = Fitness + NDSort([Population.objs,Population.cons],inf);
            end
            %% Population initialization
            P   = UniformPoint(Problem.N,Problem.D,'Latin');
            Dec = P.*repmat(Problem.upper-Problem.lower,Problem.N,1) + repmat(Problem.lower,Problem.N,1);
            Mask = double(UniformPoint(Problem.N,Problem.D,'Latin') > 0.5);
            Population = Problem.Evaluation(Dec.*Mask);
            K           = 1;  % Number of subpopulations
            Masks       = cell(1,K);
            Decs        = cell(1,K);
            Populations = cell(1,K);
            GV          = cell(1,K);
            FrontNo     = cell(1,K);
            CrowdDis    = cell(1,K);
            leader      = cell(1,K);
            index       = randperm(floor(Problem.N/K)*K);
            temp        = reshape(index,K,floor(Problem.N/K));
            for i = 1 : K
                Populations{i} = Population(temp(i,:));
                Masks{i}       = Mask(temp(i,:),:);
                Decs{i}        = Dec(temp(i,:),:);
                [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},length(Populations{i}));
                GV{i}          = UpdateGV(zeros(1,Problem.D),Masks{i},FrontNo{i});
            end
            %% Optimization
            while Algorithm.NotTerminated(Population)
                [~,rank] = sort(SubPopRank(Populations));
                for i = 1 : K
                    GV{rank(i)} = UpdateGV(GV{rank(i)},Masks{rank(i)},FrontNo{rank(i)});
                    Mating      = TournamentSelection(2,2*length(Populations{rank(i)}),FrontNo{rank(i)},-CrowdDis{rank(i)});
                    % Adaptive Variation
                    fgv = GV{rank(i)};
                    fgv(fgv<1e-2) = 0;
                    idx = kmeans(fgv',2);
                    idx = idx';
                    s1 = sum(GV{rank(i)}(idx==1))/sum(idx==1);
                    s2 = sum(GV{rank(i)}(idx==2))/sum(idx==2);
                    s = max(s1,s2);
                    delta= Problem.FE/Problem.maxFE;
                    if s < 0.5 || K==1
                        [OffDec,OffMask] = Operator1(Problem,Decs{rank(i)}(Mating,:),Masks{rank(i)}(Mating,:),Fitness);
                    else
                        leader{rank(i)}  = Leader(Masks{rank(i)},FrontNo{rank(i)},GV{rank(i)},Problem.D);
                        [OffDec,OffMask] = Operator2(Problem,Decs{rank(i)}(Mating,:),Masks{rank(i)}(Mating,:),GV{rank(i)},delta);
                    end
                    Offspring = Problem.Evaluation(OffDec.*OffMask);
                    Populations{rank(i)} = [Populations{rank(i)},Offspring];
                    Decs{rank(i)}        = [Decs{rank(i)};OffDec];
                    Masks{rank(i)}       = [Masks{rank(i)};OffMask];

                    %  Improve Environmental Selection
                    if K == 1
                        [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},CrowdDis{rank(i)}] = EnvironmentalSelection(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},floor(length(Populations{rank(i)})/2));
                    else
                        r = zeros(1,Problem.D);
                        for j = 1:K
                            if j~=i
                                fgv = GV{rank(j)};
                                fgv(fgv < 1e-2) = 0;
                                idx = kmeans(fgv',2);
                                idx = idx';
                                s1 = sum(fgv(idx==1))/sum(idx==1);
                                s2 = sum(fgv(idx==2))/sum(idx==2);
                                s = max(s1,s2);
                                if s > 0.5
                                    if s1 > s2
                                        r(idx==1) = r(idx==1)+ fgv(idx==1);
                                    else
                                        r(idx==2) = r(idx==2)+ fgv(idx==2);
                                    end
                                end
                            end
                        end
                        r(r>0)=1;
                        dis = sum(repmat(r,length(Populations{rank(i)}),1)&Masks{rank(i)},2);
                        [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},CrowdDis{rank(i)}] = EnvironmentalSelection(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},floor(length(Populations{rank(i)})/2),dis);
                    end
                end

                %% hybrid hierarchical clustering
                if mod(ceil(Problem.FE/Problem.N),10)==0
                    k = K;
                    %divide operation(Line103-118)
                    while 1
                        for i = 1 : K
                            [row,col,divisionFlag] = DivisionFlag(Masks{i},FrontNo{i});
                            if divisionFlag == 1
                                k = k+1;
                                [Populations{i},Masks{i},Decs{i},FrontNo{i},GV{i},Populations{k},Masks{k},Decs{k},FrontNo{k},GV{k}] = Divide(Populations{i},Masks{i},Decs{i},FrontNo{i},row,col);
                                leader{i}= Masks{i}(1,:);
                                leader{k}= Masks{k}(1,:);
                            end
                        end
                        K = k;
                        if divisionFlag == 0
                            break;
                        end
                    end
                    % merge operation（Line120-143）
                    fgv = GV;
                    for i = 1: K
                        fgv{i}(fgv{i} < 0.5) = 0 ;
                        fgv{i}(fgv{i} > 0) = 1;
                        mask = Masks{i}(FrontNo{i}==1,:);
                        dis = pdist2(mask,fgv{i},'hamming');
                        index = find(dis==min(dis),1);
                        fgv{i} = mask(index,:);
                    end
                    [ss,index]   = SubPopSimility(Populations,leader);
                    while K>2 && ss>0.5
                        i = index(1);
                        j = index(2);
                        [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection([Populations{i},Populations{j}],[Decs{i};Decs{j}],[Masks{i};Masks{j}],length(Populations{i}));
                        Populations(j) = [];
                        Decs(j)        = [];
                        Masks(j)       = [];
                        FrontNo(j)     = [];
                        GV(j) = [];
                        fgv(j) = [];
                        K = K-1;
                        [ss,index]   = SubPopSimility(Populations(1:K),Masks);
                    end
                    % fill operation
                    popsize = floor(Problem.N/K); %average size of each Subpopulation
                    [~,rank] = sort(SubPopRank(Populations));
                    for i = 1:K
                        len = popsize - length(Populations{rank(i)});
                        if len > 0
                            P   = UniformPoint(len,Problem.D,'Latin');
                            dec  = P.*unifrnd(repmat(Problem.lower,len,1),repmat(Problem.upper,len,1));
                            mask  = zeros(len,Problem.D);
                            F = zeros(1,Problem.D);
                            for j = 1: K
                                if j ~=i
                                    F = F + GV{rank(j)};
                                end
                            end
                            for j = 1 : len
                                mask(j,TournamentSelection(2,floor(rand*Problem.D),F)) = 1;
                            end
                            Populations{rank(i)} = [Populations{rank(i)},Problem.Evaluation(dec.*mask)];
                            Masks{rank(i)}       = [Masks{rank(i)};mask];
                            Decs{rank(i)}        = [Decs{rank(i)};dec];
                        end
                        [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},CrowdDis{rank(i)}] = EnvironmentalSelection(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},popsize);
                    end
                end
                Population = [Populations{:}];
            end
        end
    end
end