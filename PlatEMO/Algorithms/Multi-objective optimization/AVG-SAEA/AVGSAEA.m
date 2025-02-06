classdef AVGSAEA< ALGORITHM
% <2024> <multi> <real/integer> <expensive> <large>
% Adaptive variable grouping based surrogate-assisted evolutionary algorithm
% Numtrain --- 300 --- Number of train samples
% wmax     ---  20 --- Number of generations before updating models
% NumEsp   ---   2 --- Number of subPopulations
% mu       ---   5 --- Number of re-evaluated solutions in each iteration
 
%------------------------------- Reference --------------------------------
% Y. Li, X. Feng, and H. Yu. Solving high-dimensional expensive
% multiobjective optimization problems by adaptive decision variable
% grouping. IEEE Transactions on Evolutionary Computation, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yingwei Li

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [Numtrain,wmax,NumEsp,mu] = Algorithm.ParameterSet(300,20,5,5);
            
            %% Initalize the population by Latin hypercube sampling
            PopDec     = UniformPoint(Numtrain,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,Numtrain,1).*PopDec+repmat(Problem.lower,Numtrain,1));
            Arc        = Population;
            subDec     = cell(1,NumEsp);
            trainDecs  = Population.decs;
            trainObjs  =  Population.objs;
            nor_Index  = [];
            flag       = 1;
            psset      = cell(NumEsp,Problem.M);
            qsset      = cell(NumEsp,Problem.M);
            modelset   = cell(NumEsp,Problem.M);
            PopObj     = cell(1,NumEsp);
            PopCon     = cell(1,NumEsp);
            muDec      = cell(1,NumEsp);

            %% Optimization
            while Algorithm.NotTerminated(Arc)
                % Adaptive variable grouping
                [better_Cpop,bad_Cpop] = DimSelect(Arc,floor(length(Arc)/4));
                CbetterDec             = better_Cpop.decs; 
                CbadDec                = bad_Cpop.decs;
                mean_betterCDec        = mean(CbetterDec);
                mean_badCDec           = mean(CbadDec);
                difference             = abs(mean_betterCDec - mean_badCDec); 
                [~, diff_idx]          = sort(difference, 'descend'); 
                varsPerGroup           = floor(Problem.D/NumEsp);
                for i = 1 : NumEsp-1
                    nor_Index = [nor_Index,ones(1,varsPerGroup).*i];
                end
                nor_Index = [nor_Index, ones(1,Problem.D-size(nor_Index,2)).*NumEsp];
                for k = 1 : Problem.D
                    change_Index(diff_idx(k)) = nor_Index(k);
                end
                for i = 1 : NumEsp
                    trainsubDec{i} = trainDecs(:,change_Index==i);
                end
                
                %% Train models
                trainLabel = trainObjs;
                 for i = 1 : NumEsp
                     for m = 1 : Problem.M
                        [Input,ps]    = mapminmax(trainsubDec{i}',0,1);
                        Input         = Input';
                        [Output,qs]   = mapminmax(trainLabel(:,m)',0,1);
                        Output        = Output';
                        dmodel        = newrbe(Input',Output', 1 );
                        modelset{i,m} = dmodel;
                        psset{i,m}    = ps;
                        qsset{i,m}    = qs;
                    end
                 end
                [Pop,~] = EnvironmentalSelection(Arc,Problem.N);
                PopDec  = Pop.decs;
                for i = 1 : NumEsp
                     subDec{i} = PopDec(:,change_Index==i);
                     PopObj{i} = Pop.objs;
                end

                %% Evolution of subpopulations
                w = 1;
                if flag == 2 
                    % Diversity-oriented environmental selection
                    while w <= wmax
                        drawnow('limitrate');
                        for j = 1 : NumEsp
                            PopCon{j}  = calCon(PopObj{j});
                            MatingPool = TournamentSelection(2,ceil(Problem.N/2)*2,PopCon{j});
                            OffDec     = OperatorG(Problem,subDec{j}(MatingPool,:),change_Index,j);  
                            N          = size(OffDec,1);
                            OffObj     = zeros(N,Problem.M);
                                for m = 1 : Problem.M
                                    normDec       = mapminmax('apply',OffDec',psset{j,m});
                                    [normPre,~,~] = sim(modelset{j,m},normDec);
                                    OffObj(:,m)   = (mapminmax('reverse',normPre,qsset{j,m}))';
                                end
                            [PopObj{j},subDec{j}] = DVsetEnvironmentalSelection([PopObj{j};OffObj],[subDec{j};OffDec],N);
                        end
                        w = w + 1;
                    end

                else
                    % Convergence-oriented environmental selection
                    while w <= wmax
                        drawnow('limitrate');
                        for i = 1 : NumEsp
                            PopCon{i}  = calCon(PopObj{i});
                            MatingPool = TournamentSelection(2,ceil(Problem.N/2)*2,PopCon{i});
                            OffDec     = OperatorG(Problem,subDec{i}(MatingPool,:),change_Index,i);  
                            N          = size(OffDec,1);
                            OffObj     = zeros(N,Problem.M);
                            for m = 1 : Problem.M
                                normDec       = mapminmax('apply',OffDec',psset{i,m});
                                [normPre,~,~] = sim(modelset{i,m},normDec);
                                OffObj(:,m)   = (mapminmax('reverse',normPre,qsset{i,m}))';
                            end
                            allCon  = calCon([PopObj{i};OffObj]);
                            Con     = allCon(1:N);
                            newCon  = allCon(N+1:end); 
                            updated = Con > newCon;
                            subDec{i}(updated,:) = OffDec(updated,:);
                            PopObj{i}(updated,:) = OffObj(updated,:);
                        end
                        w = w + 1;
                    end
                end

                %% Merge variables in each group into complete solution
                if flag == 1
                    for i = 1 : NumEsp   
                           betterIndex = SubEnvironmentalSelection(PopObj{i},mu);
                           muDec{i}    = subDec{i}(betterIndex,:);
                           NewDec(:,change_Index==i) = muDec{i};
                    end
                elseif flag == 2
                    for j = 1 : NumEsp   
                           [~,muDec{j}] = DSelectNew(PopObj{j},subDec{j},mu);
                           NewDec(:,change_Index==j) = muDec{j};
                    end
                    flag = 1;
                else
                    finalMergeDec = zeros(Problem.N,Problem.D); 
                    for i = 1 : NumEsp
                        finalMergeDec(:,change_Index==i) = subDec{i};
                    end
                    N = size(PopObj{1},1);
                    stdSamples = zeros(N,Problem.M);
                    for m = 1 : Problem.M
                        Samples = [];
                        for i = 1 : NumEsp   
                            Samples = [Samples,PopObj{i}(:,m)];
                        end
                        stdSamples(:,m) = std(Samples,0,2);
                    end
                    meanstd     = mean(stdSamples,2);
                    [~,std_Idx] = sort(meanstd,"descend");
                    NewDec      = finalMergeDec(std_Idx(1:mu),:);
                    flag        = 1;
                end

                %% Select mu new samples
                New = Problem.Evaluation(NewDec);
                lastPopulation = Arc;
                Arc = [Arc,New]; 
                if flag == 0
                    flag = 1;
                else
                    flag = CalFlag(Arc,lastPopulation);
                end

                %% Select train data for the next iteration
                [trainDecs, trainObjs] = SelectTrainData(Arc, Numtrain, length(New));
            end
        end
    end
end

function Con = calCon(PopObj)
% Calculate the convergence of each solution

    FrontNo = NDSort(PopObj,inf);
    Con     = sum(PopObj,2);
    Con     = FrontNo'*(max(Con)-min(Con)) + Con;
end