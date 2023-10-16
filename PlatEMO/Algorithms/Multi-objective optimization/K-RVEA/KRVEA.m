classdef KRVEA < ALGORITHM
% <multi/many> <real/integer> <expensive>
% Surrogate-assisted RVEA
% alpha ---  2 --- The parameter controlling the rate of change of penalty
% wmax  --- 20 --- Number of generations before updating Kriging models
% mu    ---  5 --- Number of re-evaluated solutions at each generation

%------------------------------- Reference --------------------------------
% T. Chugh, Y. Jin, K. Miettinen, J. Hakanen, and K. Sindhya, A surrogate-
% assisted reference vector guided evolutionary algorithm for
% computationally expensive many-objective optimization, IEEE Transactions
% on Evolutionary Computation, 2018, 22(1): 129-142.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [alpha,wmax,mu] = Algorithm.ParameterSet(2,20,5);

            %% Generate the reference points and population
            [V0,Problem.N] = UniformPoint(Problem.N,Problem.M);
            V     = V0;
            NI    = 11*Problem.D-1;
            P     = UniformPoint(NI,Problem.D,'Latin');
            A2    = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            A1    = A2;  
            THETA = 5.*ones(Problem.M,Problem.D);
            Model = cell(1,Problem.M);

            %% Optimization
            while Algorithm.NotTerminated(A2)
                % Refresh the model and generate promising solutions
                A1Dec = A1.decs;
                A1Obj = A1.objs;
                for i = 1 : Problem.M
                    % The parameter 'regpoly1' refers to one-order polynomial
                    % function, and 'regpoly0' refers to constant function. The
                    % former function has better fitting performance but lower
                    % efficiency than the latter one
                    dmodel     = dacefit(A1Dec,A1Obj(:,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
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
                    index  = KEnvironmentalSelection(PopObj,V,(w/wmax)^alpha);
                    PopDec = PopDec(index,:);
                    PopObj = PopObj(index,:);
                    % Adapt referece vectors
                    if ~mod(w,ceil(wmax*0.1))
                        V(1:Problem.N,:) = V0.*repmat(max(PopObj,[],1)-min(PopObj,[],1),size(V0,1),1);
                    end
                    w = w + 1; 
                end

                % Select mu solutions for re-evaluation
                [NumVf,~] = NoActive(A1Obj,V0);
                PopNew    = KrigingSelect(PopDec,PopObj,MSE(index,:),V,V0,NumVf,0.05*Problem.N,mu,(w/wmax)^alpha); 
                New       = Problem.Evaluation(PopNew);
                A2        = [A2,New];
                A1        = UpdataArchive(A1,New,V,mu,NI); 
            end
        end
    end
end