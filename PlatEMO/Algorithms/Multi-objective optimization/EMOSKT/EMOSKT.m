classdef EMOSKT < ALGORITHM
% <2024> <multi> <real/binary> <large/none> <constrained/none> <sparse> <multitask>
% Evolutionary multi-objective optimization with sparsity knowledge transfer

%------------------------------- Reference --------------------------------
% C. Wu, Y. Tian, L. Zhang, X. Xiang, and X. Zhang. A sparsity knowledge
% transfer-based evolutionary algorithm for large-scale multitasking multi-
% objective optimization. IEEE Transactions on Evolutionary Computation,
% 2024.
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
            %% Population initialization
            TaskNum    = size(Problem.SubM,2);
            EachN      = ceil(Problem.N/TaskNum);
            [~, Maxid] = max(Problem.SubD);
            Buquan     = zeros(1,TaskNum);
            for i = 1 : TaskNum
                Buquan(i) = Problem.SubD(Maxid)-Problem.SubD(i);
            end
            % Calculate the fitness of each decision variable
            REAL    = all(Problem.encoding~=4);
            TDec    = cell(1+4*REAL,TaskNum);
            TMask   = cell(1+4*REAL,TaskNum);
            TempPop = cell(1+4*REAL,TaskNum);
            Fitness = cell(1,TaskNum);
            for j = 1:TaskNum
                Fitness{j} = zeros(1,Problem.SubD(j));
            end
            FitnessDec = cell(1,TaskNum);
            for j = 1 : TaskNum
                FitnessDec{j} = zeros(1+4*REAL,Problem.SubD(j));
            end
            for i = 1 : 1+4*REAL
                for j = 1 : TaskNum
                    if REAL
                        Dec{i,j} = unifrnd(repmat(Problem.lower(1:Problem.SubD(j))+(Problem.upper(1:Problem.SubD(j))-Problem.lower(1:Problem.SubD(j)))*((i-1)/(1+4*REAL)),Problem.SubD(j),1),...
                                           repmat(Problem.lower(1:Problem.SubD(j))+(Problem.upper(1:Problem.SubD(j))-Problem.lower(1:Problem.SubD(j)))*((i)/(1+4*REAL)),Problem.SubD(j),1));
                    else
                        Dec{i,j} = ones(Problem.SubD(j),Problem.SubD(j));
                    end
                    Mask{i,j}     = eye(Problem.SubD(j));
                    Skill{i,j}    = j*ones(Problem.SubD(j),1);
                    Solution{i,j} = [Dec{i,j}.*Mask{i,j},zeros(size(Dec{i,j},1),Buquan(j)),Skill{i,j}];
                    Initpop{i,j}  = Problem.Evaluation(Solution{i,j});
                    TDec{i,j}     = [TDec{i,j};Dec{i,j}];
                    TMask{i,j}    = [TMask{i,j};Mask{i,j}];
                    TempPop{i,j}  = [TempPop{i,j},Initpop{i,j}];
                    Fitness{j}    = Fitness{j} + NDSort([Initpop{i,j}.objs,Initpop{i,j}.cons],inf);
                    FitnessDec{j}(i,:) = NDSort([Initpop{i,j}.objs,Initpop{i,j}.cons],inf);
                end
            end
            FitnessPop = TempPop;
            for j = 1 : TaskNum
                for i = 1 : Problem.SubD(j)
                    pDecIndex{j,i} = find(FitnessDec{j}(:,i)==min(FitnessDec{j}(:,i)));
                end
            end
            % Generate initial population
            Dec  = cell(1,TaskNum);
            Mask = cell(1,TaskNum);
            Pop  = {};
            SubPopulation = {};
            FrontNo       = {};
            CrowdDis      = {};
            Skill         = {};
            Solution      = {};
            for i = 1 : TaskNum
                if REAL
                    for n = 1 : EachN
                        for d = 1 : Problem.SubD(i)
                            pDecRandIndex = pDecIndex{i,d}(randi(size(pDecIndex{i,d},1)));
                            Dec{i}(n,d)   = unifrnd(Problem.lower(d)+(Problem.upper(d)-Problem.lower(d))*((pDecRandIndex-1)/(1+4*REAL)),...
                                                    Problem.lower(d)+(Problem.upper(d)-Problem.lower(d))*((pDecRandIndex)/(1+4*REAL)));
                        end
                    end
                else
                    Dec{i} = ones(EachN,Problem.SubD(i));
                end
                Mask{i} = zeros(EachN,Problem.SubD(i));
                SamMask   = Mask{i}(1:5,:);
                [~,rank1] = sort(Fitness{i});
                index     = round([0.1,0.2,0.3,0.4,0.5]*Problem.SubD(i));
                for j = 1 : 5
                    SamMask(j,rank1(1:index(j))) = 1;
                end
                Mask{i}(1:5,:) = SamMask;
                for j = 6 : EachN
                    Mask{i}(j,TournamentSelection(2,ceil(rand*Problem.SubD(i)),Fitness{i})) = 1;
                end             
                Skill{i}    = i*ones(EachN,1);
                Solution{i} = [Dec{i}.*Mask{i},zeros(size(Dec{i},1),Buquan(i)),Skill{i}];
                Pop{i}      = Problem.Evaluation(Solution{i});
            end
            % Generate initthetamid
            initthetamid = zeros(1,TaskNum);
            for i = 1 : TaskNum
                [~,~,~,TFrontNo,~] = SparseEA_ESnouni([Pop{i},[TempPop{:,i}]],[Dec{i};vertcat(TDec{:,i})],[Mask{i};vertcat(TMask{:,i})],length([Pop{i},[TempPop{:,i}]]));
                [SubPopulation{i},Dec{i},Mask{i},FrontNo{i},CrowdDis{i}] = SparseEA_EnvironmentalSelection([Pop{i},[TempPop{:,i}]],[Dec{i};vertcat(TDec{:,i})],[Mask{i};vertcat(TMask{:,i})],EachN);
                if size(find(TFrontNo(1:5)==1),2)>0
                    initthetamid(i) = mean(index(TFrontNo(1:5)==1));
                else
                    Theta = sum(Mask{i}(FrontNo{i} ==1,:),2)';
                    initthetamid(i) = mean(Theta);
                end
            end

            %% Optimization
            NumTransUp        = floor(EachN/10)*ones(1,TaskNum);
            ALLTthetamid      = [];
            ALLTthetamid(1,:) = initthetamid;
            [SourceId,TF1]    = SourceTaskrand(Problem,SubPopulation,Dec,Mask,FrontNo,CrowdDis,EachN,Fitness,FitnessDec,pDecIndex,Buquan);
            while Algorithm.NotTerminated([SubPopulation{:}])
                for i = 1 : TaskNum
                    [SubPopulation{i},Dec{i},Mask{i},FrontNo{i},CrowdDis{i}] = OP_SparseEA(i,EachN,Problem,FrontNo{i},CrowdDis{i},Dec{i},Mask{i},Fitness{i},SubPopulation{i},Buquan(i),REAL);
                end
                [SubPopulation,Dec,Mask,FrontNo,CrowdDis,NumTransUp,ALLTthetamid] = Op_Transnodec(Problem,SubPopulation,Dec,Mask,FrontNo,CrowdDis,EachN,TF1,FitnessDec,pDecIndex,FitnessPop,SourceId,Buquan,NumTransUp,ALLTthetamid);
                [SourceId,TF1] = SourceTaskrand(Problem,SubPopulation,Dec,Mask,FrontNo,CrowdDis,EachN,Fitness,FitnessDec,pDecIndex,Buquan);
            end
        end
    end
end