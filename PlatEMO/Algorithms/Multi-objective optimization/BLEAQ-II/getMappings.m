function [psiMapping,phiMapping,lies] = getMappings(Problem,indv,archive,option)
% Construct Psi-mapping or Phi-mapping or both

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    psiApprox    = []; 
    phiApprox    = [];
    constructPsi = 0; 
    constructPhi = 0;
    
    if nargin < 4
        option = 3;
    end
    if option == 1
        constructPsi = 1;
    elseif option == 2
        constructPhi = 1;
    elseif option == 3
        constructPsi = 1;
        constructPhi = 1;
    end
    
    % Minimum number of Lower level points required for quadratic approximation of both mappings
    numberLowerLevelPointsTrain = (Problem.DU+1)*(Problem.DU+2)/2+2*(Problem.DU);
    numberLowerLevelPointsEval  = Problem.DU;
    minNumberLowerLevelPoints   = numberLowerLevelPointsTrain + numberLowerLevelPointsEval;

    PopDec = archive.decs;
    quadArchive.upper = PopDec(:,1:Problem.DU);
    quadArchive.lower = PopDec(:,Problem.DU+1:end);
    
    % Only the members close to the offspring are used for quadratic approximation
    distances = zeros(size(quadArchive.upper,1));
    for k = 1 : length(quadArchive)
        distances(k) = sum((indv - quadArchive.upper(k,:)).^2);
    end
    [~,I]  = sort(distances);
    I      = I(1:minNumberLowerLevelPoints);           
    sizeI  = length(I);
    permut = randperm(sizeI);
    J      = permut(1:sizeI-numberLowerLevelPointsEval);
    quadApproxMembers = I(J);
    % I stores the members being considered from the archive
    % quadApproxMembers stores the members from I used for a quadratic
    setDiff = setdiff(I,quadApproxMembers);
    if constructPsi == 1
        psiApprox = cell(1,Problem.DL); sumMSE = 0;
        for j = 1 : Problem.DL 
            psiApprox{j} = QuadApprox(quadArchive.lower(quadApproxMembers,j),quadArchive.upper(quadApproxMembers,:));
            %Calculating sum of MSE for all lower level variable approximations
            sumMSE = sumMSE + psiApprox{j}.mseNorm;
        end
        predictedLowerLevelVariables = zeros(length(setDiff),Problem.DL);
        for k = 1 : length(setDiff)
            for j = 1 : Problem.DL
                predictedLowerLevelVariables(k,j) = psiApprox{j}.constant + quadArchive.upper(setDiff(k),:)*psiApprox{j}.linear + quadArchive.upper(setDiff(k),:)*psiApprox{j}.sqmatrix*quadArchive.upper(setDiff(k),:)';
            end
        end
        psiMapping.function = psiApprox;
        psiMapping.sumMSE   = sumMSE;
        psiMapping.validMSE = mean(mean((predictedLowerLevelVariables - quadArchive.lower(setDiff,:)).^2));
    end
    
    if constructPhi == 1 
        PopObj        = archive.objs;
        quadArchive.llFunctionValue = PopObj(:,2);
        llFunctionDim = size(quadArchive.llFunctionValue,2);
        sumMSE        = 0;
        phiApprox     = cell(1,llFunctionDim);
        for j = 1 : llFunctionDim
            phiApprox{j} = QuadApprox(quadArchive.llFunctionValue(quadApproxMembers,j), quadArchive.upper(quadApproxMembers,:));
            sumMSE       = sumMSE+phiApprox{j}.mseNorm;
        end
        predictedLowerLevelObjective = zeros(length(setDiff),llFunctionDim);  
        for k = 1 : length(setDiff)
            for j = 1 : llFunctionDim
                predictedLowerLevelObjective(k,j) = phiApprox{j}.constant + quadArchive.upper(setDiff(k),:)*phiApprox{j}.linear + quadArchive.upper(setDiff(k),:)*phiApprox{j}.sqmatrix*quadArchive.upper(setDiff(k),:)';
            end
        end
        phiMapping.function = phiApprox;
        phiMapping.sumMSE   = sumMSE;
        phiMapping.validMSE = mean(mean((predictedLowerLevelObjective - quadArchive.llFunctionValue(setDiff,:)).^2));
    end
    
    %Checks if the offspring lies in between the training points and not outside
    lies = 1;
    for j = 1 : size(quadArchive.upper,2)
        if indv(j)>=max(quadArchive.upper(quadApproxMembers,j)) || indv(j)<=min(quadArchive.upper(quadApproxMembers,j))
            lies = 0;
        end
    end
end