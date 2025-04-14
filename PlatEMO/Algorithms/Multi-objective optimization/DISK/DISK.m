classdef DISK < ALGORITHM
% <2024> <multi/many> <real/integer> <expensive>
% Distribution-based Kriging-assisted evolutionary algorithm
% wmax  --- 60 --- Generations of evolutionary search
% alpha ---  5 --- Number of selected candidates

%------------------------------- Reference --------------------------------
% Z. Zhang, Y. Wang, G. Sun, and T. Pang. A distribution information based 
% Kriging-assisted evolutionary algorithm for expensive many-objective 
% optimization problems. IEEE Transactions on Evolutionary Computation, 
% 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: z.zhang0@csu.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            global NI mu K 
            [wmax,alpha] = Algorithm.ParameterSet(60,5);
            
            %% Initialization
            NI    = Problem.N;
            OP    = UniformPoint(NI,Problem.D,'Latin');
            A2    = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*OP+repmat(Problem.lower,NI,1));
            A1    = A2;
            THETA = 5.*ones(Problem.M,Problem.D);
            Model = cell(1,Problem.M);
              
            while Algorithm.NotTerminated(A2)
                %% Surrogate Construction
                [Model,THETA] = model_train(A2,Model,THETA);
                
                %% Learning Distribution
                [F,~]  = NDSort(A2.objs,inf);
                PopDec = A2(F==1).decs;
                if size(PopDec,1) <= 1
                    PopDec = [PopDec; A2(F==2).decs];
                end
                mu = mean(PopDec,1);
                K  = (PopDec-mu)'*(PopDec-mu)/(size(PopDec,1)-1);
                  
                %% Evolutionary Search
                OP = optimizaiton(A1,wmax,Model,Problem);
                
                %% Candidate Selection
                C = NewSelect(OP,A2,alpha,Problem);
                
                %% Adaptive Exploration
                % Judgement
                flag = 0;
                if ~isempty(C)
                    flag = judgeLS(C,A2);
                    A2   = [A2,C];
                end
                if flag == 1
                    % Surrogate Construction
                    [Model,THETA] = model_train(A2,Model,THETA);
                    % Exploration
                    [W,ideal]     = IdentifyW(A2,Problem.N,Problem.M);
                    A2            = LocalSearch(OP,W,ideal,wmax,Model,A2,Problem);
                end
                
                %% Population Update
                index = EnvironmentalSelection(A2.objs,NI);
                A1    = A2(index);
            end
        end
    end
end

function [Model,THETA] = model_train(A2,Model,THETA)
    Dec     = A2.decs;
    Obj     = A2.objs;
    Len_dec = size(Dec,2);
    Len_obj = size(Obj,2);
    for i = 1 : Len_obj
        [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
        [~,distinct2] = unique(round(Obj(:,i)*1e100)/1e100,'rows');
        distinct      = intersect(distinct1,distinct2);

        dmodel     = dacefit(Dec(distinct,:),Obj(distinct,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model{i}   = dmodel;
        THETA(i,:) = dmodel.theta;
    end
end

function [OffObj,Off_ObjMSE] = model_predict(Model,OffDec)
    N          = size(OffDec,1);
    Len_obj    = length(Model);
    OffObj     = zeros(N,Len_obj);
    Off_ObjMSE = zeros(N,Len_obj);

    for i = 1 : N
        for j = 1 : Len_obj
            [OffObj(i,j),~,Off_ObjMSE(i,j)] = predictor(OffDec(i,:),Model{j});
        end
    end
    OffObj     = real(OffObj);
    Off_ObjMSE = abs(real(Off_ObjMSE));
end

function P = optimizaiton(Population,wmax,Model,Problem)
    w      = 1;
    [N,~]  = size(Population.decs);
    P.decs = Population.decs;
    while w <= wmax    
        OffDec = OperatorGA(Problem,P.decs);
        P.decs = [P.decs;OffDec];
        [P.objs,P.objmse] = model_predict(Model,P.decs);
                
        index = SEnvironmentalSelection(P,N);  
        
        P.decs   = P.decs(index,:);
        P.objs   = P.objs(index,:);
        P.objmse = P.objmse(index,:);
             
        w = w + 1;
    end
end

function flag = judgeLS(C,A2)
    [F1,~] = NDSort(C.objs,1);
    AObj   = C(F1==1).objs;
    [F2,~] = NDSort(A2.objs,1);
    A2Obj  = A2(F2==1).objs;
    N1     = size(AObj,1);
    N2     = size(A2Obj,1);

    dominate = zeros(N1,N2);
    for i = 1 : N1
        for j = 1 : N2
            if all(AObj(i,:) <= A2Obj(j,:)) && ~all(AObj(i,:) == A2Obj(j,:))
                dominate(i,j) = true;
            end
        end
    end
    if any(any(dominate))
        flag = 0;
    else
        flag = 1;
    end
end