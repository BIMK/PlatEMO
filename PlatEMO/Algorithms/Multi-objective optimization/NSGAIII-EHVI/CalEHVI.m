function EI = CalEHVI(RealFirstObj,Obj,MSE,S,S_S)
% EHVI calculation using improtance sampling

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    nSample    = size(S,1);
    NRealFront = size(RealFirstObj,1);
    
    PopObj   = [RealFirstObj;Obj];
    [NPop,M] = size(Obj);
    Zmin     = min(PopObj,[],1);
    Zmax     = max(PopObj,[],1);
    
    % Normalization
    a   = Zmax - Zmin;
    Obj = Obj - repmat(Zmin,size(Obj,1),1);
    RealFirstObj = RealFirstObj - repmat(Zmin,size(RealFirstObj,1),1);
    Obj = Obj./repmat(a,NPop,1);
    RealFirstObj = RealFirstObj ./repmat(a,NRealFront,1);
    
    PopMSE = MSE./repmat((a.^2),size(MSE,1),1);
    sigma  = PopMSE.^0.5;
    R_S    = zeros(NRealFront,nSample);
    for i = 1 : NRealFront
        x        = sum(repmat(RealFirstObj(i,:),nSample,1)-S<=0,2) == M;
        R_S(i,x) = 1;
    end
    index     = (sum(R_S,1) == 0);
    NonDomS   = S(index,:);
    nNonDomS  = size(NonDomS,1);
    NonDomS_S = S_S(index,index);

    % Measurement of improvement
    I  = sum(NonDomS_S,2);
    EI = ones(NPop,1);
    for i = 1 : NPop
        if any(sigma(i,:)==0)
            EI(i) = 0;
        else
            P = ones(nNonDomS,1);
            for j = 1 : M
                p = normpdf(NonDomS(:,j) ,Obj(i,j),sigma(i,j));
                P = P.*p;
            end
            EI(i) = sum(I.*P);
        end
    end
end