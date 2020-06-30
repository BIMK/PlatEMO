function [OffDec,OffMask] = Operator(ParentDec,ParentMask,rbm,dae,Site,allZero,allOne,Lower,Upper,REAL)
% The operator of MOEA/DSR
%
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Parent1Mask = ParentMask(1:end/2,:);
    Parent2Mask = ParentMask(end/2+1:end,:);
    Parent1Dec  = ParentDec(1:end/2,:);
    Parent2Dec  = ParentDec(end/2+1:end,:);
    
    %% Binary variation
    if any(Site)
        other   = ~allZero & ~allOne;
        OffTemp = BinaryVariation(rbm.reduce(Parent1Mask(Site,other)),rbm.reduce(Parent2Mask(Site,other)));
        OffTemp = rbm.recover(OffTemp);
        OffMask = false(size(OffTemp,1),size(Parent1Mask,2));
        OffMask(:,other)  = OffTemp;
        OffMask(:,allOne) = true;
    else
        OffMask = [];
    end
    OffMask = [OffMask;BinaryVariation(Parent1Mask(~Site,:),Parent2Mask(~Site,:))];
    
    %% Real variation
    if REAL
        if any(Site)
            OffDec = RealVariation(dae.reduce(Parent1Dec(Site,:)),dae.reduce(Parent2Dec(Site,:)));
            OffDec = dae.recover(OffDec);
        else
            OffDec = [];
        end
        OffDec = [OffDec;RealVariation(Parent1Dec(~Site,:),Parent2Dec(~Site,:),Lower,Upper)];
    else
        OffDec = ones(size(OffMask));
    end
end

function Offspring = BinaryVariation(Parent1,Parent2)
% One point crossover and bitwise mutation

    [proC,proM] = deal(1,1);
    [N,D] = size(Parent1);
    k = repmat(1:D,N,1) > repmat(randi(D,N,1),1,D);
    k(repmat(rand(N,1)>proC,1,D)) = false;
    Offspring1    = Parent1;
    Offspring2    = Parent2;
    Offspring1(k) = Parent2(k);
    Offspring2(k) = Parent1(k);
    Offspring     = [Offspring1;Offspring2];
    Site = rand(2*N,D) < proM/D;
    Offspring(Site) = ~Offspring(Site);
end

function Offspring = RealVariation(Parent1,Parent2,Lower,Upper)
% Simulated binary crossover and polynomial mutation

    [proC,disC,proM,disM] = deal(1,20,1,20);
    [N,D] = size(Parent1);
    beta  = zeros(N,D);
    mu    = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    if nargin > 2
        Lower = repmat(Lower,2*N,1);
        Upper = repmat(Upper,2*N,1);
    else
        Lower = zeros(2*N,D);
        Upper = ones(2*N,D);
    end
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring       = min(max(Offspring,Lower),Upper);
end