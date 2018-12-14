function [etaC,Pn,allEgamma]  = OperatorAdaption(etaC,Pn,allEgamma,Egamma,nmov)
% Adapt the parameters of operators according to MRDL

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Add the current MRDL to the MRDL array
    if isempty(allEgamma)
        if ~isnan(Egamma)
            allEgamma = Egamma;
        end
        return;
    else
        if ~isnan(Egamma)
            allEgamma = [allEgamma,Egamma];
        else
            allEgamma = [allEgamma,allEgamma(end)];
        end
    end

    %% Calculate the moving average of MRDL of each generation
    MAEgamma = zeros(size(allEgamma));
    for i = 1 : length(allEgamma)
        MAEgamma(i) = mean(allEgamma(max(1,i-nmov+1):i));
    end
    
    %% Calculate the predicted current MRDL
    Y      = log(MAEgamma(1:end-1))';
    Phi    = [ones(1,length(Y));1:length(Y)]';
    Lambda = Phi\Y;
    predictEgamma = exp(Lambda(1)+Lambda(2)*length(allEgamma));
    
    %% Update the parameters of operators
    if allEgamma(end) > predictEgamma
        etaC = max(etaC-2,2);
        Pn   = 0.5*(allEgamma(end)-predictEgamma);
    else
        etaC = etaC + 2;
        Pn   = 0;
    end
end