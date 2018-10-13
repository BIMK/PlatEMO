function OffspringDec = Operator(ParentDec,lower,upper)
% Simulated binary crossover for WOF

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    ParentDec   = ParentDec([1:end,1:ceil(end/2)*2-end],:);
    [N,D]       = size(ParentDec);

    %% Simulated binary crossover
    Parent1Dec = ParentDec(1:N/2,:);
    Parent2Dec = ParentDec(N/2+1:end,:);
    beta = zeros(N/2,D);
    mu   = rand(N/2,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(20+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(20+1));
    beta = beta.*(-1).^randi([0,1],N/2,D);
    beta(rand(N/2,D)<0.5) = 1;
    beta(repmat(rand(N/2,1)>0.9,1,D)) = 1;
    OffspringDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
                    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];
    OffspringDec = max(min(OffspringDec,upper),lower);
end