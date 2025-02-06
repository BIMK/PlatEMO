classdef NSGAIIIEHVI < ALGORITHM
% <2023> <multi/many> <real> <expensive>
% NSGA-III with expected hypervolume improvement
% wmax    ---    15 --- Number of generations before updating Kriging models
% LB      ---  -0.5 --- Low bound for EHVI calculation
% UB      ---   1.2 --- Upper bound for EHVI calculation
% nSample --- 10000 --- number of samples in importance sampling for EHVI
% randp   ---   0.3 --- The parameter in generating a random number

%------------------------------- Reference --------------------------------
% Y. Pang, Y. Wang, S. Zhang, X. Lai, W. Sun, and X. Song. An expensive 
% many-objective optimization algorithm based on efficient expected 
% hypervolume improvement. IEEE Transactions on Evolutionary 
% Computation, 2023, 27(6): 1822-1836.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yong Pang (email: pangy@mail.dlut.edu.cn)

    methods
        function main(Algorithm,Problem)
            % Parameter setting
            [wmax,LB,UB,nSample, randp] = Algorithm.ParameterSet(15,-0.5,1.2,10000,0.3);
            
            % Initialization of NSGAIII
            NI            = Problem.N;
            P             = UniformPoint(NI,Problem.D,'Latin');
            Population    = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Zmin          = min(Population.objs,[],1); 
            A             = Population;
            theta         = 10.*ones(Problem.M,Problem.D);
            Model         = cell(1,Problem.M);

            % Sample and calculate hypervolume improvement
            Lowb = LB .*ones(1,Problem.M);
            Upb  = UB .*ones(1,Problem.M);
            S    = UniformPoint(nSample,Problem.M,'Latin');
            S    = S.*repmat(Upb-Lowb,nSample,1)+repmat(Lowb,nSample,1);
            S_S  = zeros(nSample,nSample);
            for i = 1 : nSample
                y        = sum(repmat(S(i,:),nSample,1)-S<=0,2) == Problem.M;  
                S_S(i,y) = 1 ;
            end 
            
            % Main loop           
            while Algorithm.NotTerminated(A)
                Dec = Population.decs;
                Obj = Population.objs;
                MSE = zeros(size(Dec,1),Problem.M);
                
                % Train kriging models
                train_X = A.decs;
                train_Y = A.objs;
                [~,distinct] = unique(round(train_X*1e6)/1e6,'rows');  
                train_X      = train_X(distinct,:);
                train_Y      = train_Y(distinct,:);
                for i = 1 : Problem.M 
                    dmodel     = dacefit(train_X,train_Y(:,i),'regpoly0','corrgauss',theta(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    theta(i,:) = dmodel.theta;
                end
                
                % Obtain current real PF
                [RealFrontNo,~] = NDSort(train_Y,1);
                RealFirstObj    = train_Y((RealFrontNo==1),:);        
                w = 1;
                while w <= wmax
                    w = w + 1;
                    OffspringDec = OperatorGA(Problem,Dec(randi(end,1,NI),:));
                    N = size(OffspringDec,1);
                    OffspringObj = zeros(N,Problem.M);
                    OffspringMSE = zeros(N,Problem.M);
                    for i = 1 : N
                        for j = 1 : Problem.M
                            [OffspringObj(i,j),~,OffspringMSE(i,j)] = predictor(OffspringDec(i,:),Model{j});
                        end
                    end
                    Zmin    = min([Zmin; OffspringObj],[],1);
                    all_Obj = [Obj;OffspringObj];
                    all_MSE = [MSE;OffspringMSE];
                    all_Dec = [Dec;OffspringDec];
                    
                    % Nodominated sorting considering the uncertainty           
                    Choose = SelectionMSE(all_Obj,all_MSE,RealFirstObj,NI,Z,Zmin);                                        
                    Dec    = all_Dec(Choose,:);
                    Obj    = all_Obj(Choose,:);
                    MSE    = all_MSE(Choose,:);
                end
     
                % Efficent EHVI calcualtion using importance sampling       
                EHVI = CalEHVI(RealFirstObj,Obj,MSE,S,S_S);
                [~,sortIndex] = sort(EHVI,'descend');

                % Diversity maintain 
                ChoseMax  = max([round(size(EHVI,1)*unifrnd(0,randp)),1]);
                FnewIndex = sortIndex(1:ChoseMax);   
                   
                % Calcualte distance
                DA     = RealFirstObj;
                dist_D = zeros(size(FnewIndex,1),size(DA,1));
                for i = 1 : size(FnewIndex,1 )
                    for j = 1 : size(DA,1)
                        dist_D(i,j) = norm(Obj(FnewIndex(i,1),:)-DA(j,:),2);
                    end
                end
                
                % Diversity Indicator
                DI = min(dist_D,[],2);  
                [ ~,SnewIndex] =max(DI);
                newIndex = [FnewIndex(SnewIndex,1)];
                PnewDec  = Dec(newIndex,:);
                PnewDec  = unique(PnewDec,'rows');               
                New  = Problem.Evaluation(PnewDec);
                A    = [A,New];
                A2   = [Population,New];
                Zmin = min(A2.objs,[],1);
                Population = EnvironmentalSelection([Population,New],Problem.N,Z,Zmin);  
            end      
        end
    end
end