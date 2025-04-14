classdef DBEMTO < ALGORITHM
% <2023> <multi> <real/integer/label/binary/permutation> <constrained>
% Double-balanced evolutionary multi-task optimization

%------------------------------- Reference --------------------------------
% K. Qiao, J. Liang, K. Yu, M. Wang, B. Qu, C. Yue, and Y. Gou. A
% self-adaptive evolutionary multi-task based constrained multi-objective
% evolutionary algorithm. IEEE Transactions on Emerging Topics in
% Computational Intelligence, 2023, 7(4): 1098-1112.
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
            Population{1} = Problem.Initialization();
            Population{2} = Problem.Initialization();

            Fm  = [0.6,0.8,1.0];
            CRm = [0.1,0.2,1.0];
            lu  = [Problem.lower;Problem.upper];
            strategy_num         = 3;
            arrayGbestChange     = ones(2,strategy_num);
            arrayGbestChangeRate = zeros(2,strategy_num);
            indexBestLN(1:2)     = 2;
            cnt = 0;
            leastSelectionPro = 0.05;
            mixPopSize        = Problem.N;
            M = Problem.M;
            numViaLN    = zeros(2,strategy_num);
            consumedFES = zeros(2,strategy_num);

            Fitness{1} = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{2} = CalFitness(Population{2}.objs);

            k = 2;

            %% Optimization
            while Algorithm.NotTerminated(Population{1})
                cnt = cnt +1;
                for m = 1 :2
                    pop{m,1} = Population{1,m}.decs;
                    pop{m,2}(:,1:M) = Population{1,m}.objs;
                end
                for i = 1 : k
                    vi = 0;
                    ui = 0;
                    arrayGbestChangeRate(i,1) = arrayGbestChange(i,1);
                    arrayGbestChangeRate(i,2) = arrayGbestChange(i,2);
                    arrayGbestChangeRate(i,3) = arrayGbestChange(i,3);

                    [~,indexBestLN(i)] = max(arrayGbestChangeRate(i,:));
                    if sum(arrayGbestChangeRate(i,:) == arrayGbestChangeRate(i,1)) == strategy_num
                        indexBestLN(i) = 1;
                    end
                    arrayGbestChange(i,:)     = 0.1 * ones(1,strategy_num);
                    arrayGbestChangeRate(i,:) = zeros(1,strategy_num);
                    consumedFES(i,:)          = ones(1,strategy_num);

                    permutation = randperm(mixPopSize);
                    if indexBestLN(i) == 1
                        arrayThird{i}  = permutation(1:ceil(leastSelectionPro*mixPopSize));
                        arraySecond{i} = permutation(ceil(leastSelectionPro*mixPopSize+1): ceil(2*leastSelectionPro*mixPopSize));
                        arrayFirst{i}  = permutation(ceil(2*leastSelectionPro*mixPopSize+1):end);
                        numViaLN(i,1)  = numViaLN(i,1) + 1;
                    elseif indexBestLN(i) == 2
                        arrayThird{i}  = permutation(1:ceil(leastSelectionPro*mixPopSize));
                        arrayFirst{i}  = permutation(ceil(leastSelectionPro*mixPopSize+1): ceil(2*leastSelectionPro*mixPopSize));
                        arraySecond{i} = permutation(ceil(2*leastSelectionPro*mixPopSize+1):end);
                        numViaLN(i,2)  = numViaLN(i,2) + 1;
                    elseif indexBestLN(i) == 3
                        arrayFirst{i}  = permutation(1:ceil(leastSelectionPro*mixPopSize));
                        arraySecond{i} = permutation(ceil(leastSelectionPro*mixPopSize+1): ceil(2*leastSelectionPro*mixPopSize));
                        arrayThird{i}  = permutation(ceil(2*leastSelectionPro*mixPopSize+1):end);
                        numViaLN(i,3)  = numViaLN(i,3) + 1;
                    end

                    rateViaLN{i}(cnt,:) = numViaLN(i,:)/sum(numViaLN(i,:));
                    consumedFES(i,:)    = consumedFES(i,:) + [length(arrayFirst{i}),length(arraySecond{i}),length(arrayThird{i})];
                    r0 = permutation;

                    %% =================Intra-task convergence strategy=================================
                    if ~isempty(arrayFirst{i})
                        if i == 1
                            MatingPool1 = TournamentSelection(2,2*length(arrayFirst{i}),Fitness{i});
                        else
                            MatingPool1 = TournamentSelection(2,2*length(arrayFirst{i}),Fitness{i});
                        end
                        valOffspring{i}(1,arrayFirst{i}) = OperatorGAhalf(Problem,Population{i}(MatingPool1));
                    end
                    %%  ====================Intra-task diversity strategy===========================
                    if ~isempty(arraySecond{i})
                        [r1,r2,r3] = gnR1R2R3(mixPopSize,  r0);
                        pop2       = pop{i,1}(arraySecond{i},:);
                        popsize2   = length(arraySecond{i});

                        index = randi([1,length(Fm)],popsize2,1);
                        F2    = Fm(index);
                        F2    = F2';
                        index = randi([1,length(CRm)],popsize2,1);
                        CR2   = CRm(index);
                        CR2   = CR2';

                        vi = pop{i,1}(r1(arraySecond{i}),:)  + F2(:, ones(1, Problem.D)) .* (pop{i,1}(r2(arraySecond{i}),:) - pop{i,1}(r3(arraySecond{i}),:));
                        vi = boundConstraint(vi, pop2, lu);

                        mask        = rand(popsize2, Problem.D) > CR2(:, ones(1, Problem.D));
                        rows        = (1 : popsize2)';
                        cols        = floor(rand(popsize2, 1) * Problem.D)+1;
                        jrand       = sub2ind([popsize2 Problem.D], rows, cols);
                        mask(jrand) = false;
                        u2          = vi;
                        u2(mask)    = pop2(mask);
                        valOffspring{i}(1,arraySecond{i}) = Problem.Evaluation(u2);
                    end
                    %% ====================Inter-task transfer strategy====================================
                    if ~isempty(arrayThird{i})
                        rand_perm   = randperm(mixPopSize);
                        MatingPool1 = rand_perm(1:length(arrayThird{i}));
                        valOffspring{i}(1,arrayThird{i}) = Population{2/i}(MatingPool1);
                    end
                end
                %% ====================Environmental selection====================================
                for i = 1 : k
                    if indexBestLN(i) ~= 3
                        if i == 1
                            [Population{i},~,Fitness{1}]    = EnvironmentalSelection([Population{i},valOffspring{2/i}],mixPopSize,i);
                            [Population{i},succ,Fitness{1}] = EnvironmentalSelection([Population{i},valOffspring{i}],mixPopSize,i);
                        elseif i == 2
                            [Population{i},~,Fitness{2}]    = EnvironmentalSelection([Population{i},valOffspring{2/i}],mixPopSize,i);
                            [Population{i},succ,Fitness{2}] = EnvironmentalSelection([Population{i},valOffspring{i}],mixPopSize,i);
                        end
                    else
                        if i == 1
                            [Population{i},succ,Fitness{1}] = EnvironmentalSelection([Population{i},valOffspring{i}],mixPopSize,i);
                        elseif i == 2
                            [Population{i},succ,Fitness{2}] = EnvironmentalSelection([Population{i},valOffspring{i}],mixPopSize,i);
                        end
                    end
                    arrayGbestChange(i,1) = arrayGbestChange(i,1) + sum(succ(arrayFirst{i}))/length(arrayFirst{i});
                    arrayGbestChange(i,2) = arrayGbestChange(i,2) + sum(succ(arraySecond{i}))/length(arraySecond{i});
                    arrayGbestChange(i,3) = arrayGbestChange(i,3) + sum(succ(arrayThird{i}))/length(arrayThird{i});
                end
            end
        end
    end
end