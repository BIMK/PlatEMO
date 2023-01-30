function [offspringsLowerLevelVariables,sumMSE,validMSE] = getLowerLevelVariableFromMapping(offsprings,psiMapping,phiMapping,Problem,archive)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % minimum size required for phi-mapping based lower level optimal variable retrieval
    minPhiSize = (Problem.DU+Problem.DL+1)*(Problem.DU+Problem.DL+2)/2 + 2*(Problem.DU) + Problem.DU;
    
    %% Determine which mapping to select in the descendant update 
    runPhi = 0;
    if (length(archive.tag1)+length(archive.tag0)) >= minPhiSize
        runPhi = ((psiMapping.sumMSE>=phiMapping.sumMSE) && (psiMapping.validMSE>=phiMapping.validMSE)) || (rand<0.25);
    end
    
	offspringsLowerLevelVariables = getLowerLevelVariableFromPsi(offsprings,psiMapping.function);
    sumMSE   = psiMapping.sumMSE;
    validMSE = psiMapping.validMSE;
        
	if runPhi == 1
		offspringsLowerLevelVariables = getLowerLevelVariableFromPhi(Problem,offsprings,offspringsLowerLevelVariables,phiMapping.function{:},[archive.tag1,archive.tag0],minPhiSize);
        sumMSE   = phiMapping.sumMSE;
        validMSE = phiMapping.validMSE;
    end
end	

function indvLowerLevelVariables = getLowerLevelVariableFromPsi(indv,psiFunction)
	DL = length(psiFunction);
	indvLowerLevelVariables = zeros(1,DL);
	for j = 1 : DL
		indvLowerLevelVariables(j) = psiFunction{j}.constant + indv*psiFunction{j}.linear + indv*psiFunction{j}.sqmatrix*indv';
	end	
end

function indvLowerLevelVariables = getLowerLevelVariableFromPhi(Problem,indv,indvLowerLevelVariables,phiFunction,archive,archiveSize)
    archiveDec       = archive.decs;
	archivePhi.upper = archiveDec(:,1:Problem.DU);
    for j = 1 : size(archivePhi.upper,1)
        distances(j) = sum((indv - archivePhi.upper(j,:)).^2);
    end
    [~,I]      = sort(distances);
    I          = I(1:archiveSize);
    archive    = archive(I);
    archiveDec = archive.decs;
    archiveObj = archive.objs;
    archiveCon = archive.cons;
    archivePhi.upper = archivePhi.upper(I,:);
    archivePhi.lower = archiveDec(:,Problem.DU+1:end);
	archivePhi.functionValue          = archiveObj(:,1);
	archivePhi.equalityConstrVals     = [];
	archivePhi.inequalityConstrVals   = archiveCon(:,1:Problem.C);	
	archivePhi.llFunctionValue        = archiveObj(:,2);
    archivePhi.llEqualityConstrVals   = [];
    archivePhi.llInequalityConstrVals = archiveCon(:,Problem.C+1:end);

    approxPhi.function = QuadApprox(archivePhi.functionValue, [archivePhi.upper archivePhi.lower]);
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
    lb      = min(archivePhi.lower);
    ub      = max(archivePhi.lower);
    options = optimset('Algorithm','sqp','Display','off');  
    [indvLowerLevelVariables,FVAL,EXITFLAG,OUTPUT] = fmincon(@(xl) -approximatedFunction(xl,indv,approxPhi.function),...
                    indvLowerLevelVariables,[],[],[],[],lb,ub,@(xl) approximatedConstraints(xl,indv,...
                    approxPhi.equalityConstr,approxPhi.inequalityConstr, approxPhi.llFunction,...
                    phiFunction, approxPhi.llEqualityConstr, approxPhi.llInequalityConstr),options);  
end

function approxFunctionValue = approximatedFunction(xl,xu,parameters)
    pop = [xu xl];
    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop'; 
end

function [c,ceq] = approximatedConstraints(llPop,ulPop,parametersEqualityConstr,parametersInequalityConstr,parametersLowerLevelFunction,parametersPhiFunction,parametersLLEqualityConstr,parametersLLInequalityConstr)
    pop = [ulPop,llPop];
    ceq = [];
    c   = [];
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
    c1       = (parametersLowerLevelFunction.constant + pop*parametersLowerLevelFunction.linear + pop*parametersLowerLevelFunction.sqmatrix*pop');
    c2       = (parametersPhiFunction.constant + ulPop*parametersPhiFunction.linear + ulPop*parametersPhiFunction.sqmatrix*ulPop');
    c(end+1) = c2 - c1;
end