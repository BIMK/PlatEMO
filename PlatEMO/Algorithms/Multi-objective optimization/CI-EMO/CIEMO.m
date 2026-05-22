classdef CIEMO < ALGORITHM
% <2026> <multi/many> <real> <expensive>
% Composite indicator-guided infilling sampling for expensive multi-objective optimization

%------------------------------- Reference --------------------------------
% H. Zhen, X. Li, W. Gong, and X. Hu. Composite indicator-guided infilling
% sampling for expensive multi-objective optimization. Swarm and
% Evolutionary Computation, 2026, 101: 102312.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

    methods
       function main(Algorithm,Problem)
            %% Population initialization
            Problem.N = 11*Problem.D-1;                                     % Initial sample size NI
            if Problem.N > 100
                Problem.N = 100;
            end
            [W, Problem.N] = UniformPoint(Problem.N,Problem.M);             % Generate reference point W
            P              = UniformPoint(Problem.N,Problem.D,'Latin');     % LHS initializes the population P, where normalized samples are generated
            Achieve        = Problem.Evaluation(repmat(Problem.upper-Problem.lower,Problem.N,1).*P+repmat(Problem.lower,Problem.N,1)); % Population evaluation, get archive Achieve
            THETA          = 5.*ones(Problem.M,Problem.D);                  % GP Parameters
            
            %% Optimization       
            while Algorithm.NotTerminated(Achieve)
                % End
                if length(Achieve) >= Problem.maxFE
                    break;
                end

                % Train GP models
                DATA  = Achieve;
                TSDec = DATA.decs;
                TSObj = DATA.objs; 
                maxnumData = 1000;                                          % 11*Problem.D-1+5;
                numTS      = size(TSDec,1);
                if size(TSDec,1) >= maxnumData
                   trainX = TSDec(numTS-maxnumData+1:end,:);
                   trainY = TSObj(numTS-maxnumData+1:end,:);
                else
                   trainX = TSDec;
                   trainY = TSObj;
                end
                [trainX, trainY] = dsmerge(trainX, trainY);                 % Combining repetitive training points
                for i = 1 : Problem.M
                    dmodel     = dacefit(trainX,trainY(:,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    THETA(i,:) = dmodel.theta;
                end 

                % NSGAIII search
                [PopDec,PopObj,MSE] = NSGAIII_opt(Achieve,Model,W,Problem);

                % Candidate Seletction
                PopNew = CS(PopDec,PopObj,Achieve.decs,Achieve.objs);
                
                % Expensive evaluation
                New = [];
                if ~isempty(PopNew)                                         % Avoiding empty samples and removing duplicate points
                    [~,ib] = intersect(PopNew,Achieve.decs,'rows');
                    PopNew(ib,:) = [];
                    if ~isempty(PopNew)
                        New = Problem.Evaluation(PopNew);                   % Evaluation of sampling points obtained
                    else
                        [~,ib]       = intersect(PopDec,Achieve.decs,'rows');
                        PopDec(ib,:) = [];
                        [a1, a2]     = size(PopDec);
                        index        = randi(a2);
                        PopNew       = PopDec(index,:);
                        New          = Problem.Evaluation(PopNew);  
                    end
                end
                Achieve = cat(2,Achieve,New);                               % Adding evaluation points to the archive
            end
       end
    end
end