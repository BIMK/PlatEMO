function EIMEGO(Global)
% <algorithm> <A-G>
% Expected Improvement Matrix-based Infill Criteria for Expensive
% Multiobjective Optimization
% InfillCriterionIndex --- 1 --- infill criterion index number

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Dawei Zhan

    %% Parameter setting
    % 1 for Euclidean distance-based EIM criterion (default)
    % 2 for Maximin distance-based EIM criterion
    % 3 for Hypervolume-based EIM criterion
    InfillCriterionIndex = Global.ParameterSet(1);

    %% Generate the initial design points
    % number of design variables
    D = Global.D;
    % number of objective functions
    M = Global.M;
    % number of initial design points
    N  = 11*D-1;
    % generate initial design points using Latin Hypercube sampling
    PopDec = repmat(Global.upper-Global.lower,N,1).*lhsamp(N,D)+repmat(Global.lower,N,1);
    % calculate initial design points
    Population   = INDIVIDUAL(PopDec);

    %% Optimization
    while Global.NotTermination(Population)
        % scale the objective values to 0 and 1
        PopDec = Population.decs;
        PopObj = Population.objs;
        N  = size(PopDec,1);       
        PopObjScaled = (PopObj-repmat(min(PopObj),N,1))./(repmat(max(PopObj),N,1)-repmat(min(PopObj),N,1));   
        % bulid Kriging models for all the objective functions
        KrigingModel = cell(1,M);
        for i = 1 : M
            KrigingModel{i}= dacefit(PopDec,PopObjScaled(:,i),'regpoly0','corrgauss',1*ones(1,D),0.001*ones(1,D),1000*ones(1,D));
        end       
        % select one candidate with the maximum EIM value using GA
        PopDec     = InfillSamplingEIM(Global,KrigingModel,PopObjScaled,InfillCriterionIndex);
        Population = [Population,INDIVIDUAL(PopDec)];
    end
end