function [rbm,dae,allZero,allOne] = ModelTraining(Mask,Dec,REAL,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Basic variable setup
    allZero = all(~Mask,1);
    allOne  = all(Mask,1);
    other   = ~allZero & ~allOne;
    
    %% Non-dominated sorting for evaluation
    [FrontNo,MaxFront] = NDSort(Dec.*Mask,ones(size(Dec,1),1));
    
    %% Variable evaluation
    validMask = Mask(:,other);
    validDec  = Dec(:,other);
    
    if ~isempty(validMask)
        % Calculate variable contributions
        contribution = abs(validDec).*validMask;
        meanContrib  = mean(contribution,1);
        
        % Normalize meanContrib to [0, 1]
        meanContribNorm = (meanContrib - min(meanContrib)) / (max(meanContrib) - min(meanContrib) + eps);
        
        % Calculate variable ranking scores
        [~,rankIdx] = sort(meanContrib,'descend');
        rankScore   = zeros(1,length(meanContrib));
        rankScore(rankIdx) = linspace(1,0,length(rankIdx));
        
        % Analyze non-dominated solutions
        bestMask = validMask(FrontNo==1,:);
        if ~isempty(bestMask)
            successRate = mean(bestMask,1);
        else
            successRate = mean(validMask,1);
        end
        
        % Calculate comprehensive scores with normalized meanContrib
        varScore = 0.4*rankScore + 0.4*successRate + 0.2*meanContribNorm;
        
        % Adaptive threshold
        progress  = Problem.FE/Problem.maxFE;
        threshold = mean(varScore) * (1 + 0.2*(1-progress));
        Type      = varScore > threshold;
    else
        Type     = [];
        varScore = [];
    end
    
    %% Determine hidden layer size
    nValid = sum(other);
    if nValid == 0
        K = 1;
    else
        % Base K on valid variables and progress
        baseK = min(sum(Type), round(sqrt(nValid)));
        K     = min(max(round(baseK * (1-0.3*progress)), 1), size(Mask,1));
    end
    
    %% Train networks
    if nValid > 0
        rbm = RBM(nValid,K,10,1,0,0.5,0.1);
        rbm.train(Mask(:,other));
    else
        rbm = [];
    end
    
    if REAL
        dae = DAE(size(Dec,2),K,10,size(Dec,1),0.5,0.5,0.1);
        dae.train(Dec);
    else
        dae = [];
    end
end