function [eliteIndivLS,llEliteIndivLS,localSearch] = DoLocalSearch(Problem,initialIndivLS,archive,StoppingCriteria,gen)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Flag.localSearch = 0; 
    Flag.recheck     = 0; 
    Flag.update      = 0; 
    Flag.replace     = 0;
    Flag.dontDoPhi   = 0;
    
    numberLowerLevelPointsTrain = (Problem.DU+1)*(Problem.DU+2)/2+2*(Problem.DU);
    numberLowerLevelPointsEval  = Problem.DU;
    minSizePsiLS = numberLowerLevelPointsTrain + numberLowerLevelPointsEval;
    minSizePhiLS = (Problem.DU+Problem.DL+1)*(Problem.DU+Problem.DL+2)/2 + 2*(Problem.DU) + Problem.DU;

    terminationFlag   = (Problem.FE > Problem.maxFE) || (StoppingCriteria == 1);
    frequency         = Problem.N;
    frequencyOffsetLS = frequency - 1;
    frequencyOffsetRC = floor(frequency/2 - 1);
    
    if length(archive.tag1) >= minSizePsiLS
        if terminationFlag == 1
            Flag.localSearch = 1; % Exact local search
        elseif mod(gen+frequencyOffsetLS,frequency) == 0
            Flag.localSearch = 2; % Default local search
        end
    end

    IndivDec     = initialIndivLS.decs;
    eliteIndiv   = IndivDec(:,1:Problem.DU);
    llEliteIndiv = IndivDec(:,Problem.DU+1:end);
    ulDimMax     = Problem.upper(:,1:Problem.DU);
    ulDimMin     = Problem.lower(:,1:Problem.DU);
    llDimMax     = Problem.upper(:,Problem.DU+1:end);
    llDimMin     = Problem.lower(:,Problem.DU+1:end);
    range        = 0.05;
    
    doExactLSFlag = Flag.localSearch;
    if length([archive.tag1,archive.tag0]) < minSizePhiLS
        Flag.dontDoPhi = 1;
    end
    
    % Obtain psi and phi mappings around eliteIndiv
    [psiMapping,phiMapping,~]   = getMappings(Problem,eliteIndiv,archive.tag1);
    functionLowerLevelVariables = psiMapping.function;
    functionLowerLevelObjective = phiMapping.function{:};
      
    localSearch.psiMSE = psiMapping.validMSE; localSearch.phiMSE = phiMapping.validMSE;
    % Adaptive selection Which mapping to pick during execution
    if (psiMapping.validMSE <= phiMapping.validMSE) || Flag.dontDoPhi == 1 || (rand <=0.5)% ÓÃpsi-mapping
        if doExactLSFlag == 1       % Run an exact local search
            localSearchMethod = 1;	% psi w/ EXACT obj. func.
        else% Run an approximate local search (default)
            localSearchMethod = 2;	% psi w/ Approx. obj. func.
        end
    else % phi-mapping
        if doExactLSFlag == 1
            localSearchMethod = 3;	% phi w/ EXACT obj. func.
        else
            localSearchMethod = 4;	% phi w/ Approx. obj. func.
        end       
    end
    
    % Data preparation
    archiveLS    = [archive.tag1,archive.tag0];
    archiveLSDec = archiveLS.decs;
    upper        = archiveLSDec(:,1:Problem.DU);
    for j = 1 : size(upper,1)
        distances(j) = sum((eliteIndiv - upper(j,:)).^2);
    end
    [~,I] = sort(distances);
    
    if localSearchMethod == 2
        archiveConsidered = (Problem.DU+1)*(Problem.DU+2)/2+2*(Problem.DU)+Problem.DU;
        % If run psi w/ approx. obj. func. 
        % need to approximate F(x) and G(x) with Tag 1 member
        archiveLSDec = archiveLS(I(1:archiveConsidered)).decs;
        archiveLSObj = archiveLS(I(1:archiveConsidered)).objs;
        archiveLSCon = archiveLS(I(1:archiveConsidered)).cons;
        archivePsi.upper = upper(I(1:archiveConsidered),:);
        archivePsi.lower = archiveLSDec(:,Problem.DU+1:end);
        archivePsi.functionValue        = archiveLSObj(:,1);
        archivePsi.equalityConstrVals   = [];
        archivePsi.inequalityConstrVals = archiveLSCon(:,1:Problem.C);
        approxPsi.function              = QuadApprox(archivePsi.functionValue, [archivePsi.upper]);
        if size(archivePsi.equalityConstrVals,2) ~= 0
	        for i = 1 : size(archivePsi.equalityConstrVals,2)  
                approxPsi.equalityConstr{i} = QuadApprox(archivePsi.equalityConstrVals(:,i), [archivePsi.upper]);
	        end
	    else
	        approxPsi.equalityConstr = [];
        end
	    if size(archivePsi.inequalityConstrVals,2) ~= 0
	        for i = 1 : size(archivePsi.inequalityConstrVals,2)
                approxPsi.inequalityConstr{i} = QuadApprox(archivePsi.inequalityConstrVals(:,i), [archivePsi.upper]);
	        end
	    else
	        approxPsi.inequalityConstr = [];
        end
    end
    
    if localSearchMethod == 4
        % If run phi w/ approx. obj. func. 
        % need to approximate F(x,y), G(x,y), f(x,y), g(x,y) w/ random
        % member
        archiveConsidered = (Problem.DU+Problem.DL+1)*(Problem.DU+Problem.DL+2)/2 + 2*(Problem.DU) + Problem.DU;
        archiveLSDec      = archiveLS(I(1:archiveConsidered)).decs;
        archiveLSObj      = archiveLS(I(1:archiveConsidered)).objs;
        archiveLSCon      = archiveLS(I(1:archiveConsidered)).cons;
        archivePhi.upper  = upper(I(1:archiveConsidered),:);
        archivePhi.lower  = archiveLSDec(:,Problem.DU+1:end);
        archivePhi.functionValue          = archiveLSObj(:,1);
        archivePhi.equalityConstrVals     = [];
        archivePhi.inequalityConstrVals   = archiveLSCon(:,1:Problem.C);
        archivePhi.llFunctionValue        = archiveLSObj(:,2);
        archivePhi.llEqualityConstrVals   = [];
        archivePhi.llInequalityConstrVals = archiveLSCon(:,Problem.C+1:end);
        approxPhi.function = QuadApprox(archivePhi.functionValue, [archivePhi.upper, archivePhi.lower]);        
		if ~isempty(archivePhi.equalityConstrVals)
		    for i = 1 : size(archivePhi.equalityConstrVals,2)		            
                approxPhi.equalityConstr{i} = QuadApprox(archivePhi.equalityConstrVals(:,i), [archivePhi.upper archivePhi.lower]);
		    end
		else
		    approxPhi.equalityConstr = [];
		end

		if ~isempty(archivePhi.inequalityConstrVals)
		    for i = 1 : size(archivePhi.inequalityConstrVals,2)
               approxPhi.inequalityConstr{i} = QuadApprox(archivePhi.inequalityConstrVals(:,i), [archivePhi.upper archivePhi.lower]);            
		    end
		else
		    approxPhi.inequalityConstr = [];
		end

		% Phi-function with Approximated F(x,y) and f(x,y)
        % Needs to approximate f(x,y), g(x,y) & h(x,y)
        approxPhi.llFunction = QuadApprox(archivePhi.llFunctionValue, [archivePhi.upper archivePhi.lower]);
 
        if ~isempty(archivePhi.llEqualityConstrVals)
            for i = 1 : size(archivePhi.llEqualityConstrVals,2)
                approxPhi.llEqualityConstr{i} = QuadApprox(archivePhi.llEqualityConstrVals(:,i), [archivePhi.upper archivePhi.lower]);
            end
        else
            approxPhi.llEqualityConstr = [];
        end

        if ~isempty(archivePhi.llInequalityConstrVals)
            for i = 1 : size(archivePhi.llInequalityConstrVals,2)
                approxPhi.llInequalityConstr{i} = QuadApprox(archivePhi.llInequalityConstrVals(:,i), [archivePhi.upper archivePhi.lower]);
            end
        else
            approxPhi.llInequalityConstr = [];
        end
    end

    if localSearchMethod == 1
        % Exact Psi-function approximation based local search 
        options = optimset('Algorithm','sqp','Display','off');
        [lb,ub] = createLocalSearchBound([eliteIndiv],[ulDimMin],[ulDimMax],range);
        [eliteIndivLS,~,EXITFLAG,OUTPUT] = fmincon(@(x) -approximatedFunctionPsi(x,llDimMin,llDimMax,functionLowerLevelVariables,Problem),[eliteIndiv],[],[],[],[],lb,ub,@(x) approximatedConstraintsPsi(x,llDimMin,llDimMax,functionLowerLevelVariables,Problem),options);
        llEliteIndivLS                   = llEliteIndiv;
        localSearch.method               = 'Psi';
        localSearch.termination          = EXITFLAG;
        localSearch.functionEvaluation   = OUTPUT.funcCount;
        return;
    end
    
    if localSearchMethod == 2
         options = optimset('Algorithm','sqp','Display','off');      
        % psi-mapping based local search w/ approximated obj. func. 
        lb = ulDimMin;
        ub = ulDimMax;
        [eliteIndivLS,~,EXITFLAG,OUTPUT] = fmincon(@(x) -approximatedFunction(x,approxPsi.function),[eliteIndiv],[],[],[],[],lb,ub,@(x) approximatedConstraints(x,approxPsi.equalityConstr,approxPsi.inequalityConstr),options);
        llEliteIndivLS                   = llEliteIndiv;
        localSearch.method               = 'Approx';
        localSearch.termination          = EXITFLAG;
        localSearch.functionEvaluation   = OUTPUT.funcCount;
        return;
    end
        
    if localSearchMethod == 3
    	% Exact Phi-function approximation based local search 
        options = optimset('Algorithm','sqp','Display','off');
        [lb,ub] = createLocalSearchBound([eliteIndiv llEliteIndiv],[ulDimMin llDimMin],[ulDimMax llDimMax],range);
        [eliteIndivFull,~,EXITFLAG,OUTPUT] = fmincon(@(x) -approximatedFunctionPhi(x,Problem),[eliteIndiv llEliteIndiv],[],[],[],[],lb,ub,@(x) approximatedConstraintsPhi(x,functionLowerLevelObjective,Problem),options);
        eliteIndivLS                   = eliteIndivFull(1:Problem.DU);
        llEliteIndivLS                 = eliteIndivFull(Problem.DU+1:end);
        localSearch.method             = 'Phi';
        localSearch.termination        = EXITFLAG;
        localSearch.functionEvaluation = OUTPUT.funcCount;
        return;
    end
    
    if localSearchMethod == 4
        % phi-mapping based local search w/ approximated obj. func.
        options = optimset('Algorithm','sqp','Display','off');
        lb      = [ulDimMin llDimMin];
        ub      = [ulDimMax llDimMax];
        [eliteIndivFull,~,EXITFLAG,OUTPUT] = fmincon(@(x) -approximatedFunction(x,approxPhi.function),[eliteIndiv llEliteIndiv],[],[],[],[],lb,ub,@(x) approximatedConstraintsPhi2(x,approxPhi.equalityConstr,approxPhi.inequalityConstr, approxPhi.llFunction, functionLowerLevelObjective, approxPhi.llEqualityConstr, approxPhi.llInequalityConstr, Problem.DU, Problem.DL),options);
        eliteIndivLS                       = eliteIndivFull(1:Problem.DU);
        llEliteIndivLS                     = eliteIndivFull(Problem.DU+1:end);
        localSearch.method                 = 'ApproxPhi';
        localSearch.termination            = EXITFLAG;
        localSearch.functionEvaluation     = OUTPUT.funcCount;    
        return;
    end
end

function functionValue = approximatedFunctionPsi(xu,llDimMin,llDimMax,psiFunction,Problem)
    for j = 1 : size(psiFunction,2)
        xl(j) = psiFunction{j}.constant + xu*psiFunction{j}.linear + xu*psiFunction{j}.sqmatrix*xu';
    end
    % check if predicted lower level optimal solution is outside the bound
    xl            = checkLimits(xl,llDimMin,llDimMax);
	Population    = Problem.Evaluation([xu,xl]);
    Obj           = Population.objs;
    functionValue = Obj(:,1);
end

function functionValue = approximatedFunctionPhi(pop,Problem)
    xu            = pop(:,1:Problem.DU);
    xl            = pop(:,Problem.DU+1:end);
	Population    = Problem.Evaluation([xu,xl]);
    Obj           = Population.objs;
    functionValue = Obj(:,1);
end
        
function approxFunctionValue = approximatedFunction(pop,parameters)
    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop';
end

function [c,ceq] = approximatedConstraintsPsi(xu,llDimMin,llDimMax,psiFunction,Problem)
    for j = 1 : size(psiFunction,2)
        xl(j) = psiFunction{j}.constant + xu*psiFunction{j}.linear + xu*psiFunction{j}.sqmatrix*xu';
    end
    % Check if predicted lower level optimal solution is outside the bound
    xl  = checkLimits(xl,llDimMin,llDimMax);
    Population = Problem.Evaluation([xu,xl]);
    Con = Population.cons;
    c   = Con(:,1:Problem.C);
    ceq = [];
end
    
function [c,ceq] = approximatedConstraintsPhi(pop,parametersPhiFunction,Problem)
    ulPop = pop(:,1:Problem.DU);
    llPop = pop(:,Problem.DU+1:end);
    Population = Problem.Evaluation(pop);
    c      = Population.cons;
    Obj    = Population.objs;
    c1     = Obj(:,2);
    ceq    = [];  
    n      = length(c);
    c2     = (parametersPhiFunction.constant + ulPop*parametersPhiFunction.linear + ulPop*parametersPhiFunction.sqmatrix*ulPop');
    c(n+1) = c2 - c1;
end

function [c,ceq] = approximatedConstraints(pop,parametersEqualityConstr,parametersInequalityConstr)
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

function [c,ceq] = approximatedConstraintsPhi2(pop,parametersEqualityConstr,parametersInequalityConstr,parametersLowerLevelFunction,parametersPhiFunction,parametersLLEqualityConstr,parametersLLInequalityConstr,dimULPop,dimLLPop)
    ulPop = pop(:,1:dimULPop);
    llPop = pop(:,dimULPop+1:dimULPop+dimLLPop);
    ceq   = [];
    c     = [];
    if ~isempty(parametersEqualityConstr)
        for i = 1 : length(parametersEqualityConstr)
            ceq(i) = parametersEqualityConstr{i}.constant + pop*parametersEqualityConstr{i}.linear + pop*parametersEqualityConstr{i}.sqmatrix*pop';
        end
    end
    if ~isempty(parametersLLEqualityConstr)
        for i = 1 : length(parametersLLEqualityConstr)
            ceq(end+1) = parametersLLEqualityConstr{i}.constant + pop*parametersLLEqualityConstr{i}.linear + pop*parametersLLEqualityConstr{i}.sqmatrix*pop';
        end
    end
    if ~isempty(parametersInequalityConstr)
        for i = 1 : length(parametersInequalityConstr)
            c(i) = parametersInequalityConstr{i}.constant + pop*parametersInequalityConstr{i}.linear + pop*parametersInequalityConstr{i}.sqmatrix*pop';
        end
    end
    if ~isempty(parametersLLInequalityConstr)
        for i = 1 : length(parametersLLInequalityConstr)
            c(end+1) = parametersLLInequalityConstr{i}.constant + pop*parametersLLInequalityConstr{i}.linear + pop*parametersLLInequalityConstr{i}.sqmatrix*pop';
        end
    end
    c1 = (parametersLowerLevelFunction.constant + pop*parametersLowerLevelFunction.linear + pop*parametersLowerLevelFunction.sqmatrix*pop');
    c2 = (parametersPhiFunction.constant + ulPop*parametersPhiFunction.linear + ulPop*parametersPhiFunction.sqmatrix*ulPop');
    c(end+1) = c2 - c1;
end
    
function offsprings = checkLimits(offsprings,DimMin,DimMax)
    numOffsprings = size(offsprings,1);
    dimMinMatrix  = DimMin(ones(1,numOffsprings),:);
    offsprings(offsprings<dimMinMatrix) = dimMinMatrix(offsprings<dimMinMatrix);
    dimMaxMatrix  = DimMax(ones(1,numOffsprings),:);
    offsprings(offsprings>dimMaxMatrix) = dimMaxMatrix(offsprings>dimMaxMatrix);
end

function [LB,UB] = createLocalSearchBound(Indv,DimMin,DimMax,range)
	diff = range.*(DimMax - DimMin);
	LB   = Indv - diff;
	LB   = checkLimits(LB,DimMin,DimMax);
	UB   = Indv + diff;
	UB   = checkLimits(UB,DimMin,DimMax);
end