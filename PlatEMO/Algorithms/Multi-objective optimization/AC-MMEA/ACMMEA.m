classdef ACMMEA < ALGORITHM
% <2023> <multi> <real/integer> <large/none> <multimodal> <sparse>
% Adaptive merging and coordinated offspring generation based multi-modal multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% X. Wang, T. Zheng, and Y. Jin. Adaptive merging and coordinated offspring
% generation in multi-population evolutionary multi-modal multi-objective
% optimization. Proceedings of the International Conference on Data-driven
% Optimization of Complex Systems, 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Xiangyu Wang (email: xiangyu.wang@uni-bielefeld.de)

    methods
        function main(Algorithm,Problem)
            %% Population initialization
            Dec  = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            Mask = zeros(Problem.N,Problem.D);
            GV   = ones(1,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),GV)) = 1;
                GV(Mask(i,:)==1) = GV(Mask(i,:)==1)+1;
            end
            Population  = Problem.Evaluation(Dec.*Mask);
            [slst] = Clustering(Population.decs, 20, [Problem.lower,Problem.upper], Problem.D);
            K=size(slst,2);
            Masks       = cell(1,K);
            Decs        = cell(1,K);
            Populations = cell(1,K);
            GV          = cell(1,K);
            FrontNo     = cell(1,K);
            CrowdDis    = cell(1,K);
            Fitness     = cell(1,K); 
            for i = 1 : K
                Populations{i} = Population(slst{i});
                Masks{i}       = Mask(slst{i},:);
                Decs{i}        = Dec(slst{i},:);
                [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},length(Populations{i}));
                GV{i}          = UpdateGV(zeros(1,Problem.D),Masks{i},FrontNo{i});
            end
            StageIFlag = 1;
            Timea      = 0;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Problem.FE < 0.3 * Problem.maxFE 
                    [~,rank] = sort(SubPopRank(Populations));
                    for i = 1 : K
                        GV{rank(i)}          = UpdateGV(GV{rank(i)},Masks{rank(i)},FrontNo{rank(i)});
                        Mating               = TournamentSelection(2,2*length(Populations{rank(i)}),FrontNo{rank(i)},-CrowdDis{rank(i)});
                        [OffDec,OffMask]     = Operator(Problem,Decs{rank(i)}(Mating,:),Masks{rank(i)}(Mating,:),GV{rank(i)}, StageIFlag);
                        Offspring            = Problem.Evaluation(OffDec.*OffMask);
                        Populations{rank(i)} = [Populations{rank(i)},Offspring];
                        Decs{rank(i)}        = [Decs{rank(i)};OffDec];
                        Masks{rank(i)}       = [Masks{rank(i)};OffMask];
                        if i > 1
                            for j = 1 : i-1
                                [~,fs(rank(j))] = min(mean(Populations{rank(j)}.objs,2));
                            end
                            R = zeros(1,Problem.D);
                            for j = 1 : i-1
                                R = R + Masks{rank(j)}(fs(rank(j)),:);
                            end
                            R(R>0) = 1;
                            dis = sum(repmat(R,length(Populations{rank(i)}),1)&Masks{rank(i)},2);
                            [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},CrowdDis{rank(i)}] = EnvironmentalSelection(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},floor(Problem.N/K),dis);
                        else
                            [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},CrowdDis{rank(i)}] = EnvironmentalSelection(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},floor(Problem.N/K));
                        end
                    end
                    if mod(ceil(Problem.FE/Problem.N),50)==0 
                        [Populations,Masks,Decs,GV,K]=SubPopSimility(Populations,Masks,Decs,GV);
                        FrontNo  = cell(1,K);
                        CrowdDis = cell(1,K);
                        for i = 1 : K
                            [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},floor(Problem.N/K));
                        end
                    end
                else 
                    if Timea == 0
                        for i = 1 : K
                            [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},Fitness{rank(i)}] = EnvironmentalSelectionS(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},length(Populations{rank(i)}));
                            GV{rank(i)} = UpdateGV(zeros(1,Problem.D),Masks{rank(i)},FrontNo{rank(i)});
                        end
                    end
                    Timea = Timea + 1;
                    StageIFlag = 0;
                    for i=1:K
                        GV{rank(i)}          = UpdateGV(GV{rank(i)},Masks{rank(i)},FrontNo{rank(i)});
                        Mating               = TournamentSelection(2,2*length(Populations{rank(i)}),FrontNo{rank(i)},Fitness{rank(i)});
                        [OffDec,OffMask]     = Operator(Problem,Decs{rank(i)}(Mating,:),Masks{rank(i)}(Mating,:),GV{rank(i)}, StageIFlag);
                        Offspring            = Problem.Evaluation(OffDec.*OffMask);
                        Populations{rank(i)} = [Populations{rank(i)},Offspring];
                        Decs{rank(i)}        = [Decs{rank(i)};OffDec];
                        Masks{rank(i)}       = [Masks{rank(i)};OffMask]; 
                        [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},Fitness{rank(i)}] = EnvironmentalSelectionS(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},floor(Problem.N/K));

                    end
                    if mod(ceil(Problem.FE/Problem.N),50)==0
                        [Populations,Masks,Decs,GV,K]=SubPopSimility(Populations,Masks,Decs,GV);
                        FrontNo  = cell(1,K);
                        CrowdDis = cell(1,K);
                        for i = 1 : K
                            [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},floor(Problem.N/K));
                        end
                        Timea    = 0;
                        [~,rank] = sort(SubPopRank(Populations));
                    end
                end
                Population = [Populations{:}];
            end
        end
    end
end