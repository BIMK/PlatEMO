classdef PIMD < ALGORITHM
% <2024> <multi/many> <real> <expensive>
% Probability and mapping crowding distance
% wmax --- 15 --- Maximum repeat time of offspring generation
% eta  ---  5 --- Maximum number of samples

%------------------------------- Reference --------------------------------
% Y. Li, W. Li, Y. Zhao, and S. Li. An infill sampling criterion based on
% improvement of probability and mapping crowding distance for expensive
% multi/many-objective optimization. Engineering Applications of Artificial
% Intelligence, 2024, 133: 108616.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yang Li (email: liyangwust@163.com)
    
    methods
        function main(Algorithm,Problem)
           %% Parameter setting
            [wmax,eta] = Algorithm.ParameterSet(15, 5);

           %% Initialize population
            NI         = Problem.N;
            P          = UniformPoint(NI,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            [W,~] = UniformPoint(Problem.N,Problem.M);
            Zmin  = min(Population.objs,[],1);
            Arc   = Population;
            THETA = 5.*ones(Problem.M,Problem.D);
            Model = cell(1,Problem.M);

           %% Optimization
            while Algorithm.NotTerminated(Arc)
                FrontNo = NDSort(Arc.objs,NI);
                CompPop = Arc(FrontNo<=1); % Comparison population
                Dec     = Population.decs; Obj = Population.objs; Mse = zeros(size(Obj,1),size(Obj,2));
                Lp      = Shape_Estimate(Population,Problem.N); % PF shape estimation 
                
                train_X = Arc.decs;
                train_Y = Arc.objs;
                [~,distinct] = unique(round(train_X*1e6)/1e6,'rows');
                train_X      = train_X(distinct,:);
                train_Y      = train_Y(distinct,:);
                
                for i = 1 : Problem.M % Train surrogates
                    dmodel     = dacefit(train_X,train_Y(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    THETA(i,:) = dmodel.theta;
                end
                
                % Model-based optimization
                w = 1;
                while w <= wmax
                    w = w + 1;
                    OffspringDec = OperatorGA(Problem,Dec(randi(end,1,length(Dec)),:));
                    N = size(OffspringDec,1);
                    OffspringObj = zeros(N,Problem.M);
                    OffspringMse = zeros(N,Problem.M);
                    for i = 1 : N
                        for j = 1 : Problem.M
                            [OffspringObj(i,j),~,OffspringMse(i,j)] = predictor(OffspringDec(i,:),Model{j});
                        end
                    end
                    all_Obj = [Obj;OffspringObj]; all_Dec = [Dec;OffspringDec]; all_Mse = [Mse;OffspringMse];
                    Zmin    = min([Zmin; OffspringObj],[],1);
                    Next    = EnvironmentalSelection_NSGAIII([Obj;OffspringObj],NI,W,Zmin);
                    Dec     = all_Dec(Next,:); Obj = all_Obj(Next,:); Mse = all_Mse(Next,:);
                end
                
                %Infill sampling
                deleteindex = find(sum(Mse,2)<10^-15);
                Dec(deleteindex,:) = [];
                Obj(deleteindex,:) = [];
                Mse(deleteindex,:) = [];
                Comobj = CompPop.objs;
                Pr_dominate = zeros(1,size(Obj,1));
                for i = 1 : size(Obj,1)
                    Pr = normcdf(Comobj,repmat(Obj(i,:),size(Comobj,1),1), repmat(sqrt(max(Mse(i,:),0)),size(Comobj,1),1) );
                    Pr_dominate(i) = -1*max(prod(Pr,2)); % Improvement of probability
                end
                
                [NormObj,NormCompare] = normalization(Obj,CompPop.objs);
                tran_NormObj          = Obj_tran(NormObj,Lp);
                tran_NormCompare      = Obj_tran(NormCompare,Lp);
                dist = pdist2(tran_NormCompare,tran_NormObj); 
                DI   = -min(dist,[],1); %Mapping crowding distance
                
                newObj      = [DI;Pr_dominate]';
                [FrontNo,~] = NDSort(newObj,1);
                PnewDec     = Dec((FrontNo==1),:);  % Find solutions in the first front
                [PnewDec,selectind] = unique(PnewDec,'rows');
                len = size(PnewDec,1);
                if size(PnewDec,1) > eta
                    pr      = Pr_dominate(FrontNo==1);
                    [~,ind] = min(pr(selectind));
                    temp    = randperm(len,eta);
                    if ~ismember(ind,temp)
                        temp(1) = ind;
                    end
                    PnewDec = PnewDec(temp,:);
                end
                New  = Problem.Evaluation(PnewDec); %Expensive evaluation
                Arc  = [Arc,New];
                Zmin = min(Arc.objs,[],1);
                Next = EnvironmentalSelection_NSGAIII(Arc.objs,NI,W,Zmin);
                Population = Arc(Next); % Update population
            end
        end
    end
end

function tran_Obj = Obj_tran(PopObj,Lp)
    PopObj   = PopObj+10^-6;
    tran_Obj = PopObj./(sum(abs(PopObj).^Lp,2)).^(1/Lp);  
end