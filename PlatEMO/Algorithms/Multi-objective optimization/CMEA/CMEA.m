classdef CMEA < ALGORITHM
% <2026> <multi> <real/integer> <large/none> <multimodal> <sparse>
% Cooperative multi-population evolutionary algorithm

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, Z. Gao, and X. Zhang. Exploring equivalent sparse
% optimal solutions of large-scale multi-modal multi-objective optimization
% problems using cooperative evolution. Information Sciences, 2026, 749:
% 123539.
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
            Dec  = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1)); 
            Mask = zeros(Problem.N,Problem.D);
            cols = size(Mask,2); 
            for i = 1 : Problem.N
                Mask(i,randi(cols,1)) = 1;
            end            
            Population  = Problem.Evaluation(Dec.*Mask);
            K           = 8;  
            Masks       = cell(1,K);
            Decs        = cell(1,K);
            Populations = cell(1,K);
            GV          = cell(1,K);
            FrontNo     = cell(1,K);
            CrowdDis    = cell(1,K);
            index       = randperm(floor(Problem.N/K)*K);
            temp        = reshape(index,K,floor(Problem.N/K));
            for i = 1 : K
                Populations{i} = Population(temp(i,:));
                Masks{i}       = Mask(temp(i,:),:);
                Decs{i}        = Dec(temp(i,:),:);
                [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},length(Populations{i}));
                GV{i}          = UpdateGV(zeros(1,Problem.D),Masks{i},FrontNo{i});
            end
            endingFlag = 0;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Subpopulation evolution
                [~,rank] = sort(SubPopRank(Populations));
                for i = 1 : K                
                    GV{rank(i)}        = UpdateGV(GV{rank(i)},Masks{rank(i)},FrontNo{rank(i)});
                    Mating             = TournamentSelection(2,2*length(Populations{rank(i)}),FrontNo{rank(i)},-CrowdDis{rank(i)});
                    Mating1            = Mating(1:floor(end/2)-1);
                    Mating2            = Mating(floor(end/2):end);
                    [OffDec1,OffMask1] = Operator(Problem,Decs{rank(i)}(Mating1,:),Masks{rank(i)}(Mating1,:),GV{rank(i)});
                    [OffDec2,OffMask2] = Operator2(Problem,Decs{rank(i)}(Mating2,:),Masks{rank(i)}(Mating2,:),sum(OffMask1,1)./size(OffMask1,1));
                    OffDec             = [OffDec1; OffDec2];
                    OffMask            = [OffMask1; OffMask2];

                    Offspring            = Problem.Evaluation(OffDec.*OffMask);
                    Populations{rank(i)} = [Populations{rank(i)},Offspring];
                    Decs{rank(i)}        = [Decs{rank(i)};OffDec];
                    Masks{rank(i)}       = [Masks{rank(i)};OffMask];
                    if i > 1 && Problem.FE/Problem.maxFE < 0.5
                        for j = 1 : i-1
                            [~,fs(rank(j))] = min(mean(Populations{rank(j)}.objs,2));
                        end
                        R = zeros(1,Problem.D);
                        for j = 1 : i-1
                            R = R + Masks{rank(j)}(fs(rank(j)),:);
                        end
                        R(R>0) = 1;
                        dis    = sum(repmat(R,length(Populations{rank(i)}),1)&Masks{rank(i)},2);
                        [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},CrowdDis{rank(i)}] = EnvironmentalSelection(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},floor(Problem.N/K),dis);
                    else
                        [Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},FrontNo{rank(i)},CrowdDis{rank(i)}] = EnvironmentalSelection(Populations{rank(i)},Decs{rank(i)},Masks{rank(i)},floor(Problem.N/K));
                    end
                end
                % Merge-and-divide operation
                if mod(ceil(Problem.FE/Problem.N),50)==0 && (Problem.FE < (Problem.maxFE - ceil(Problem.N/K) - 1))
                	[~,best]     = SubPopRank(Populations);
                    divisionFlag = all(best==1);
                    [ss,index]   = SubPopSimility(Populations,Masks);
                    if ss > 0.5
                        K = K-1;
                        i = index(1);
                        j = index(2);
                        [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection([Populations{i},Populations{j}],[Decs{i};Decs{j}],[Masks{i};Masks{j}],floor(Problem.N/K));
                        Populations(j) = [];
                        Decs(j)        = [];
                        Masks(j)       = [];
                        GV(j)          = [];
                        endingFlag     = endingFlag + 1;
                    elseif divisionFlag == 1
                        for i = 1 : K
                            [Populations{i},Decs{i},Masks{i},FrontNo{i},CrowdDis{i}] = EnvironmentalSelection(Populations{i},Decs{i},Masks{i},floor(Problem.N/(K+1)));
                        end
                        K    = K + 1;
                        Dec  = unifrnd(repmat(Problem.lower,floor(Problem.N/K),1),repmat(Problem.upper,floor(Problem.N/K),1));
                        Mask = zeros(floor(Problem.N/K),Problem.D);
                        F    = zeros(1,Problem.D);
                        for i = 1: K-1
                            F = F + GV{i};
                        end
                        for i = 1 : floor(Problem.N/K)
                            Mask(i,TournamentSelection(2,floor(rand*Problem.D),F)) = 1;
                        end
                        Populations{K} = Problem.Evaluation(Dec.*Mask);
                        Masks{K}       = Mask;
                        Decs{K}        = Dec;
                        GV{K}          = zeros(1,Problem.D);
                        [Populations{K},Decs{K},Masks{K},FrontNo{K},CrowdDis{K}] = EnvironmentalSelection(Populations{K},Decs{K},Masks{K},length(Populations{K}));
                        GV{K}          = UpdateGV(zeros(1,Problem.D),Masks{K},FrontNo{K});
                    end
                end
                Population = [Populations{:}];
            end
        end
    end
end