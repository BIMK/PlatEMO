function Sigma = UpdateCMA(X,Sigma,gen)
% Update the CMA model

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    n = size(X,2);
    
    %% Calculate the CMA parameters
    mu    = 4 + floor(3*log(n));
    mu1   = floor(mu/2);
    w     = log((mu+1)/2) - log(1:mu1);
    w     = w./sum(w);
    mueff = 1./sum(w.^2);
    cs    = (mueff+2)./(n+mueff+5);
    ds    = 1 + 2*max(0,sqrt((mueff-1)./(n+1))-1) + cs;
    cc    = (4+mueff/n)./(n+4+2*mueff/n);
    c1    = 2./((n+1.3).^2+mueff);
    cmu   = min(1-c1,2*(mueff-2+1/mueff)./((n+2).^2+mueff)); % Modified
    ENI   = sqrt(n)*(1-1/4/n+1/21/n^2);
    
    %% Update the CMA model
    y           = (X(1:mu1,:)-repmat(Sigma.x,mu1,1))/Sigma.sigma;
    yw          = w*y;
    Sigma.x     = Sigma.x + Sigma.sigma*yw;
    Sigma.ps    = (1-cs)*Sigma.ps + sqrt(cs*(2-cs)*mueff)*Sigma.C^(-1/2)*yw';
    hs          = norm(Sigma.ps)./sqrt(1-(1-cs).^(2*(gen+1))) < (1.4+2/(n+1))*ENI;
    deltahs     = 1 - hs; % Modified
    Sigma.pc    = (1-cc)*Sigma.pc + hs*sqrt(cc*(2-cc)*mueff)*yw;
    Sigma.sigma = Sigma.sigma*exp(cs/ds*(norm(Sigma.ps)/ENI-1));
    Sigma.C     = (1-c1-cmu)*Sigma.C + c1*(Sigma.pc'*Sigma.pc+deltahs*Sigma.C) + cmu*y'*diag(w)*y;
    Sigma.C     = triu(Sigma.C) + triu(Sigma.C,1)'; % Enforce symmetry
    
    %% Reset the CMA model if possible
    [B,D] = eig(Sigma.C);
    diagD = diag(D);
    diagC = diag(Sigma.C);
    ConditionCov  = max(diagD) > 1e14*min(diagD);
    NoEffectCoord = any(Sigma.x==Sigma.x+0.2*Sigma.sigma*sqrt(diagC)');
    NoEffectAxis  = all(Sigma.x==Sigma.x+0.1*Sigma.sigma*sqrt(diagD(mod(gen,n)+1))*B(:,mod(gen,n)+1)');
    TolXUp        = any(Sigma.sigma*sqrt(diagC)>1e4);
    if ConditionCov || NoEffectCoord || NoEffectAxis || TolXUp
        Sigma = struct('s',[],'x',[],'sigma',0.5,'C',eye(n),'pc',0,'ps',0);
    end
end