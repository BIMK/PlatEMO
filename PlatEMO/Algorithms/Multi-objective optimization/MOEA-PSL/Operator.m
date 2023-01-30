function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,rbm,dae,Site,allZero,allOne)
% The operator of MOEA/PSL
%
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
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
        OffTemp = BinaryCrossover(rbm.reduce(Parent1Mask(Site,other)),rbm.reduce(Parent2Mask(Site,other)));
        OffTemp = rbm.recover(OffTemp);
        OffMask = false(size(OffTemp,1),size(Parent1Mask,2));
        OffMask(:,other)  = OffTemp;
        OffMask(:,allOne) = true;
    else
        OffMask = [];
    end
    OffMask = [OffMask;BinaryCrossover(Parent1Mask(~Site,:),Parent2Mask(~Site,:))];
    OffMask = BinaryMutation(OffMask);
    
    %% Real variation
    if any(Problem.encoding~=4)
        if any(Site)
            OffDec = RealCrossover(dae.reduce(Parent1Dec(Site,:)),dae.reduce(Parent2Dec(Site,:)));
            OffDec = dae.recover(OffDec);
        else
            OffDec = [];
        end
        OffDec = [OffDec;RealCrossover(Parent1Dec(~Site,:),Parent2Dec(~Site,:))];
        OffDec = RealMutation(OffDec,Problem.lower,Problem.upper);
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(size(OffMask));
    end
end

function Offspring = BinaryCrossover(Parent1,Parent2)
% Uniform crossover

    k = rand(size(Parent1)) < 0.5;
    Offspring1    = Parent1;
    Offspring2    = Parent2;
    Offspring1(k) = Parent2(k);
    Offspring2(k) = Parent1(k);
    Offspring     = [Offspring1;Offspring2];
end

function Offspring = BinaryMutation(Offspring)
% Bitwise mutation

    Site = rand(size(Offspring)) < 1/size(Offspring,2);
    Offspring(Site) = ~Offspring(Site);
end

function Offspring = RealCrossover(Parent1,Parent2)
% Simulated binary crossover

    disC  = 20;
    [N,D] = size(Parent1);
    beta  = zeros(N,D);
    mu    = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
end

function Offspring = RealMutation(Offspring,Lower,Upper)
% Polynomial mutation

    disM  = 20;
    [N,D] = size(Offspring);
    Lower = repmat(Lower,N,1);
    Upper = repmat(Upper,N,1);
    Site  = rand(N,D) < 1/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring       = min(max(Offspring,Lower),Upper);
end