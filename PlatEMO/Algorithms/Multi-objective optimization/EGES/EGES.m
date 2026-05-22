classdef EGES < ALGORITHM
% <2025> <multi/many> <real> <large> <expensive>
% Efficient grouping evolutionary search

%------------------------------- Reference --------------------------------
% H. Zhen, W. Gong, L. Wang, and X. Hu. Surrogate-assisted efficient
% grouping evolutionary search for expensive large-scale multi-objective
% optimization. IEEE Transactions on Evolutionary Computation, 2025.
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
            Problem.N = 11*Problem.D-1;                                                        % Initial sample size
            if Problem.N >100
                Problem.N = 100;
            end
            P       = UniformPoint(Problem.N,Problem.D,'Latin');                               % LHS initializes the population P, where normalized samples are generated
            Achieve = Problem.Evaluation(repmat(Problem.upper-Problem.lower,Problem.N,1).*P+repmat(Problem.lower,Problem.N,1)); % Population evaluation, get archive Achieve

            %% Optimization
            while Algorithm.NotTerminated(Achieve)
                % Termination
                if length(Achieve) >= Problem.maxFE
                    break;
                end

                % DATA
                DATA  = Achieve;
                TSDec = DATA.decs;
                TSObj = DATA.objs; 
                maxnumData = 300;                                                              
                numTS      = size(TSDec,1);
                if size(TSDec,1)>=maxnumData
                   trainX = TSDec(numTS-maxnumData+1:end,:);
                   trainY = TSObj(numTS-maxnumData+1:end,:);
                else
                   trainX = TSDec;
                   trainY = TSObj;
                end
                [trainX, trainY] = dsmerge(trainX, trainY);                                     % Merge duplicate training points

                % RBF model
                pair   = pdist2(trainX, trainX);                                                % Calculate the distance matrix
                D_max  = max(max(pair, [], 2));                                                 % Find the maximum distance
                spread = D_max * (Problem.D * Problem.N) ^ (-1 / Problem.D);                    % Calculate empirical parameters
                net    = newrbe(transpose(trainX), transpose(trainY), spread);
                Model  = @(x) sim(net,x');                                                      % Build surrogate model 
                
                % GES
                [PopDec,PopObj] = GES(Achieve,Model,Problem);
                
                % CIS
                PopNew = CIS(PopDec,PopObj,Achieve.decs,Achieve.objs);
                PopNew = unique(PopNew, 'rows');

                % Expensive evaluation
                New = [];
                if ~isempty(PopNew)                                                             % Avoid empty sampling and delete duplicate points
                    [~,ib]       = intersect(PopNew,Achieve.decs,'rows');
                    PopNew(ib,:) = [];
                    if ~isempty(PopNew)
                        New = Problem.Evaluation(PopNew);                                       % Evaluated sampling points
                    else
                        [~,ib]       = intersect(PopDec,Achieve.decs,'rows');
                        PopDec(ib,:) = [];
                        [a1, a2]     = size(PopDec);
                        index        = randi(a1);
                        PopNew       = PopDec(index,:);
                        New          = Problem.Evaluation(PopNew);
                    end
                end
                Achieve = cat(2,Achieve,New);                                                   % Update archive
            end
       end
    end
end