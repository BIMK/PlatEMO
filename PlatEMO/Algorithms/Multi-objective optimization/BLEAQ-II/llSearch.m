function [eliteIndiv,tag] = llSearch(Problem,ulPopDec,llPopmember,llmemberVariance, llPop)
% Obtain the upper member corresponds to the best lower member

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% SQP
    epsilonZero = 1e-6;
    tag         = 0;
    % Initialize the lower population
    llPopSize = (Problem.DL+1)*(Problem.DL+2)/2+Problem.DL;
    if ~isempty(llPopmember)
        for i = 1 : llPopSize - 1
            llPopDec(i,:) = Problem.lower(Problem.DU+1:end) + rand(1, Problem.DL).*(Problem.upper(Problem.DU+1:end)-Problem.lower(Problem.DU+1:end));
        end
        llPopDec(llPopSize,:) = llPopmember;
    else
        for i = 1 : llPopSize
            llPopDec(i,:) = Problem.lower(Problem.DU+1:end) + rand(1, Problem.DL).*(Problem.upper(Problem.DU+1:end)-Problem.lower(Problem.DU+1:end));
        end
        llPopmember = llPopDec(llPopSize,:);
    end
    llPopulation = Problem.EvaluationLower([repmat(ulPopDec,llPopSize,1),llPopDec]);
    llPopCon     = llPopulation.cons;
    llPopObj     = llPopulation.objs;
    if sum(sum(isnan(llPopCon)))>0 || sum(sum(isnan(llPopObj)))>0
        llPopObj = llPopObj(~isnan(llPopObj));
        llPopCon = llPopCon(~isnan(llPopCon));
        llPopCon = reshape(llPopCon,[llPopSize,length(llPopCon)/llPopSize]);
    end
    % Construct the Lower quadratic approximations of objective function and linear approximations of constraints
    approx.function       = QuadApprox(llPopObj, llPopDec);
    approx.equalityConstr = [];
    if size(llPopCon,2) ~= 0
        for i = 1 : size(llPopCon,2)
            approx.inequalityConstr{i} = QuadApprox(llPopCon(:,i), llPopDec);
        end
    else
        approx.inequalityConstr = [];
    end
    
    % Quadratic functions of linear constraints are optimized using sequence quadratic programming
    options = optimset('Display','off','TolX',1e-10,'TolFun',1e-10);
    llUpper = Problem.upper(Problem.DU+1:end);
    llLower = Problem.lower(Problem.DU+1:end);
    [eliteIndivllDec,~,~,~,~] = fmincon(@(x) approximatedFunction(x,approx.function),llPopmember,[],[],[],[],llLower,llUpper,@(x) approximatedConstraints(x,approx.equalityConstr,approx.inequalityConstr),options);
    eliteIndiv   = Problem.EvaluationLower([ulPopDec, eliteIndivllDec]);
    llPopulation = [llPopulation, eliteIndiv];
    
    eliteIndivDec               = eliteIndiv.decs;
    eliteFunctionValue          = eliteIndiv.objs;
    elitellFunctionValue        = eliteFunctionValue(~isnan(eliteFunctionValue));
    eliteInequalityConstrVals   = eliteIndiv.cons;
    elitellInequalityConstrVals = eliteInequalityConstrVals(~isnan(eliteInequalityConstrVals));
    elitellEqualityConstrVals   = [];
    
    f       = approximatedFunction(eliteIndivDec(:, Problem.DU+1:end),approx.function);
    [c,ceq] = approximatedConstraints(eliteIndivllDec,approx.equalityConstr,approx.inequalityConstr);
    d       = sqrt((f-elitellFunctionValue)^2+sum((c-elitellInequalityConstrVals).^2)+sum((ceq-elitellEqualityConstrVals).^2));
    if d < epsilonZero
        tag = 1;
        eliteIndiv = eliteIndiv.dec(Problem.DU+1:end);
        return;
    end
    if tag == 0
        llPopDec = zeros(Problem.N, Problem.DL);
        if ~isempty(llPop)
            if size(llPop,1) > Problem.N
                r = randperm(size(llPop,1));
                r = r(1:Problem.N);
                llPopDec = llPop(r,:);
            else
                for i = 1 : Problem.N
                    llPopDec(i,:) = Problem.lower(Problem.DU+1:end) + rand(1, Problem.DL).*(Problem.upper(Problem.DU+1:end)-Problem.lower(Problem.DU+1:end));
                end
                llPopDec(1:size(llPop,1),:) = llPop;
            end
        else
            for i = 1 : Problem.N
                llPopDec(i,:) = Problem.lower(Problem.DU+1:end) + rand(1, Problem.DL).*(Problem.upper(Problem.DU+1:end)-Problem.lower(Problem.DU+1:end));
            end
            llPopDec(Problem.N+1,:) = eliteIndivllDec;
        end
    end
    llPopulation = Problem.EvaluationLower([repmat(ulPopDec, size(llPopDec,1),1), llPopDec]);
    FElower      = 0;
    alpha_init   = sum(var(llPopDec));
    
    %% Optimization
    while FElower < Problem.maxFElower
        % Select parents and generate offspring
        MatingPool  = TournamentSelection(2,3,CalFitness(Problem.C,llPopulation));
        ParentDec   = llPopulation(MatingPool).decs;
        llOffDec    = OperatorPCX(ParentDec(:,Problem.DU+1:end),Problem.lower(Problem.DU+1:end),Problem.upper(Problem.DU+1:end));
        llOffspring = Problem.EvaluationLower([repmat(ulPopDec,size(llOffDec,1),1),llOffDec]);
        FElower     = FElower + length(llOffspring);
        % Select r members with better adaptability
        llPopulation = EnvironmentalSelection(Problem,llPopulation,llOffspring);
        llPopDec = llPopulation.decs;
        alpha = sum(var(llPopDec(:,Problem.DU+1:end)))/alpha_init;
        if alpha < 1e-4
            tag = 1;
            [~,best]   = min(CalFitness(Problem.C,llPopulation));
            eliteIndiv = llPopulation(best).dec(Problem.DU+1:end);
            return;
        end
    end
    [~,best]   = min(CalFitness(Problem.C,llPopulation));
    eliteIndiv = llPopulation(best).dec(Problem.DU+1:end);
    
end

function approxFunctionValue = approximatedFunction(pop, parameters)
    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop';
end

function [c, ceq] = approximatedConstraints(pop, parametersEqualityConstr, parametersInequalityConstr)
    if ~isempty(parametersEqualityConstr)
        for i = 1 : length(parametersEqualityConstr)
            ceq(i) = parametersEqualityConstr{i}.constant + pop*parametersEqualityConstr{i}.linear + pop*parametersEqualityConstr{i}.sqmatrix*pop';
        end
    else
        ceq = [];
    end

    if ~isempty(parametersInequalityConstr)
        for i = 1 : length(parametersInequalityConstr)
            c(i) = parametersInequalityConstr{i}.constant + pop*parametersInequalityConstr{i}.linear + pop*parametersInequalityConstr{i}.sqmatrix*pop';
        end
    else
        c = [];
    end
end