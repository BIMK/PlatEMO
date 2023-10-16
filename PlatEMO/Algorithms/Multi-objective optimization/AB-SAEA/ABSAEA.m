classdef ABSAEA < ALGORITHM
% <multi/many> <real/integer> <expensive>
% Adaptive Bayesian based surrogate-assisted evolutionary algorithm
% alpha ---  2 --- The parameter controlling the rate of change of penalty
% wmax  --- 20 --- Number of generations before updating Kriging models
% mu    ---  5 --- Number of re-evaluated solutions at each generation

%------------------------------- Reference --------------------------------
% X. Wang, Y. Jin, S. Schmitt S, and M. Olhofer. An adaptive Bayesian
% approach to surrogate-assisted evolutionary multi-objective optimization.
% Information Sciences, 2020, 519: 317-331.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [alpha,wmax,mu] = Algorithm.ParameterSet(2,20,5);

            %% Generate the reference points and population
            [V0,Problem.N] = UniformPoint(Problem.N,Problem.M);
            V     = V0;
            V1    = V0;
            NI    = 11*Problem.D-1;
            P     = UniformPoint(NI,Problem.D,'Latin');
            A1    = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));  
            L     = 11*Problem.D-1+25; % L = 300;
            THETA = 5.*ones(Problem.M,Problem.D);
            Model = cell(1,Problem.M);

            %% Optimization
            while Algorithm.NotTerminated(A1)
                % Refresh the model and generate promising solutions delete the multiple data
                [~,index]  = unique(A1.decs,'rows');
                A1Dec = A1.decs;
                A1Dec = A1Dec(index,:);                
                A1Obj = A1.objs;
                A1Obj = A1Obj(index,:);
                Numdata=size(A1Dec,1);
                % Limit the size of the data
                if Numdata > L
                    [FrontNo,~] = NDSort(A1Obj,Numdata);
                    [~,index] = sort(FrontNo);
                    A1Dec1 = A1Dec(index(1:floor(L/2)), :);
                    A1Obj1 = A1Obj(index(1:floor(L/2)), :);                    
                    A1Dec2 = A1Dec(index(1:L-floor(L/2)), :);
                    A1Obj2 = A1Obj(index(1:L-floor(L/2)), :);
                    index = randperm(size(A1Dec2,1));
                    A1Dec = [A1Dec1;A1Dec2(index(1:L-floor(L/2)),:)];
                    A1Obj = [A1Obj1;A1Obj2(index(1:L-floor(L/2)),:)];
                end
                for i = 1 : Problem.M
                    % The parameter 'regpoly1' refers to one-order polynomial
                    % function, and 'regpoly0' refers to constant function. The
                    % former function has better fitting performance but lower
                    % efficiency than the latter one  
                    [mS, mY]   = dsmerge(A1Dec, A1Obj(:,i));
                    dmodel     = dacefit(mS, mY,'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    THETA(i,:) = dmodel.theta;                  
                end
                PopDec = A1Dec;
                w      = 1;
                while w <= wmax
                    drawnow('limitrate');
                    OffDec = OperatorGA(Problem,PopDec);
                    PopDec = [PopDec;OffDec];
                    [N,~]  = size(PopDec);
                    PopObj = zeros(N,Problem.M);
                    MSE    = zeros(N,Problem.M);
                    for i = 1: N
                        for j = 1 : Problem.M
                            [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),Model{j});
                        end
                    end
                    Selection = FSelection(PopObj,V,(w/wmax)^alpha,2);
                    PopDec = PopDec(Selection,:);
                    PopObj = PopObj(Selection,:);
                    MSE=MSE(Selection,:);
                    % Adapt referece vectors
                    if(mod(w, ceil(wmax*0.1)) == 0)
                        V = V0.*repmat(max(A1Obj,[],1)-min(A1Obj,[],1),size(V0,1),1);
                    end
                    w = w + 1;
                end

                % Select mu solutions for re-evaluation acquisiion function
                a=-0.5*cos(Problem.FE*pi/Problem.maxFE)+0.5;
                b=0.5*cos(Problem.FE*pi/Problem.maxFE)+0.5;
                [MMSE,~]=max(MSE,[],1);
                [MPopObj,~]=max(PopObj,[],1);
                fit=PopObj./MPopObj*b+MSE./MMSE*a;
                % Select by the reference vectors
                if a>0.5
                    Flag=2;
                else
                    Flag=1;
                end                
                Selection = FSelection(fit,V1,(Problem.FE/Problem.maxFE)^alpha,Flag);
                PopNew1 = PopDec(Selection,:);
                if size(PopNew1,1)<mu
                    PopNew=PopNew1(:,1:Problem.D);
                elseif size(PopNew1,1)>=mu
                    PopNew=PopNew1(randperm(size(PopNew1,1),5),1:Problem.D);
                end
                
                % Adapt referece vectors
                if(mod(Problem.FE, ceil(Problem.maxFE*0.1)) == 0)
                    V1 = V0.*repmat(max(A1Obj,[],1)-min(A1Obj,[],1),size(V0,1),1);
                end
                PopNew1 = Problem.Evaluation(PopNew);
                A1 = [A1 PopNew1];
            end
        end
    end
end