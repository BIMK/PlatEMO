classdef EMMOEA < ALGORITHM
% <multi> <real/integer> <expensive>
% Expensive multi-/many-objective evolutionary algorithm
% gmax --- 10 --- Number of generations before updating Kriging Kmodels

%------------------------------- Reference --------------------------------
% S. Qin, C. Sun, Q. Liu, and Y. Jin, A performance indicator-based infill
% criterion for expensive multi-/many-objective optimization, IEEE
% Transactions on Evolutionary Computation, 2023, 27(4): 1085-1099.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by ShufenQin

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            gmax = Algorithm.ParameterSet(10);
            
            %% Initialization
            [V,~]  = UniformPoint(Problem.N,Problem.M);
            NI     = 100;
            Plhs   = UniformPoint(NI,Problem.D,'Latin');
            TS     = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1) .* Plhs + repmat(Problem.lower,NI,1));
            lob    = 1e-5.*ones(1,Problem.D);
            upb    = 100.*ones(1,Problem.D);
            THETA  = 5.*ones(Problem.M,Problem.D);
            THETA0 = 5.*ones(1,Problem.D);
            Kmodel = cell(1,Problem.M);
            
            %% Optimization
            PopDec = TS.decs;
            while Algorithm.NotTerminated(TS)
                %% Train a GP model for each objective
                TSDec = TS.decs;
                TSObj = TS.objs;
                for i = 1 : Problem.M
                    % The parameter 'regpoly1' refers to one-order polynomial
                    % function, and 'regpoly0' refers to constant function. The
                    % former function has better fitting performance but lower
                    % efficiency than the latter one
                    dmodel     = dacefit(TSDec,TSObj(:,i),'regpoly0','corrgauss',THETA(i,:),lob,upb);
                    Kmodel{i}  = dmodel;
                    THETA(i,:) = dmodel.theta;
                end
                
                %% Model-based optimization
                g = 1;
                while g <= gmax
                    drawnow();
                    OffDec = OperatorGA(Problem,PopDec);
                    PopDec = [PopDec;OffDec];
                    N  = size(PopDec,1);
                    PopObj = zeros(N,Problem.M);
                    MSE    = zeros(N,Problem.M);
                    for i = 1: N
                        for j = 1 : Problem.M
                            [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),Kmodel{j});
                        end
                    end
                    index  = KrigingSelection(PopObj,V);
                    PopDec = PopDec(index,:);
                    g = g + 1;
                end
                
                %% Surrogate management
                % Normalization
                TSObj = (TSObj- min(TSObj,[],1))./(max(TSObj,[],1)-min(TSObj,[],1));
                Z = min(TSObj,[],1);
                
                % Calcute the performance indicator
                dc  = pdist2(TSObj,Z,'euclidean');
                ddt = pdist2(TSObj,TSObj,'euclidean');
                ddt(ddt==0) = inf;
                dd = min(ddt,[],2);
                IP = dc-dd;
                IPmin = min(IP);
                
                % Train a GP model for the performance indicator
                IPmodel = dacefit(TSDec,IP,'regpoly0','corrgauss',THETA0,lob,upb);
                THETA0  = IPmodel.theta;
                PreIP   = zeros(size(PopDec,1),1);
                MSEIP   = zeros(size(PopDec,1),1);
                for i = 1: size(PopDec,1)
                    [PreIP(i),~,MSEIP(i)] = predictor(PopDec(i,:),IPmodel);
                end
                s     = sqrt(MSEIP);
                lamda = zeros(size(PopDec,1),1);
                EIP   = zeros(size(PopDec,1),1);
                for i = 1:size(PopDec,1)
                    lamda(i) = (IPmin-PreIP(i))/s(i);
                    EIP(i) = (IPmin-PreIP(i))*normcdf(lamda(i))+s(i)*normpdf(lamda(i));
                end
                [~,maxind] = max(EIP);
                Popreal = PopDec(maxind,:);
                
                % Delete duplicated solutions
                TSDec    = [TS.decs;Popreal];
                [~,index] = unique(TSDec,'rows');
                if length(index) == size(TSDec,1)
                    PopNew = Problem.Evaluation(Popreal);
                    TS     = [TS,PopNew];
                end
                
                % Update the non-dominated solution set
                TSnd = TS((NDSort(TS.objs,1)==1));
                
                %% The next population
                TSndObj = TSnd.objs;
                TSndDec = TSnd.decs;
                TSndObj = (TSndObj - min(TSndObj,[],1))./(max(TSndObj,[],1)-min(TSndObj,[],1));
                
                AngleND = acos(1-pdist2(TSndObj,V,'cosine'));
                [ValueAngle,associate] = min(AngleND,[],2);
                TSsec = [];
                for i = unique(associate)'
                    current     = find(associate==i);
                    [~,minindc] = min(ValueAngle(current));
                    TSsect      = TSndDec(current(minindc),:);
                    TSsec       = [TSsec;TSsect];
                end
                PopDec = [PopDec;TSsec];
                PopDec = unique(PopDec,'rows');
            end
        end
    end
end