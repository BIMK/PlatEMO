function [OffDec,OffMask] = Reproduction2(Problem,ParentDec,ParentMask,NSV,SV,REAL)
% The operator of MOEA/CKF

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Parent1Mask = ParentMask(1:end/2,:);
    Parent2Mask = ParentMask(end/2+1:end,:);
    [N,D] = size(Parent1Mask);
    
    %% Binary variation
    OffMask = false(N,D);
    OffMask(:,SV) = BinaryCrossover(Parent1Mask(:,SV),Parent2Mask(:,SV));
    OffMask(:,SV) = BinaryMutation(OffMask(:,SV));
    
    OffMask(:,NSV) = BinaryCrossover(Parent1Mask(:,NSV),Parent2Mask(:,NSV));
    OffMask(:,NSV) = BinaryMutation(OffMask(:,NSV));
    
    %% Real variation
    if REAL
        Parent1 = ParentDec(1:floor(end/2),:);
        Parent2 = ParentDec(floor(end/2)+1:floor(end/2)*2,:);             
        OffDec  = zeros(N,D);
        for i = 1 : N
            NonZero = OffMask(i,:);
            OffDec(i,NonZero)  = SBXhalf(Parent1(i,NonZero),Parent2(i,NonZero));
            OffDec(i,NonZero)  = PM(Problem,OffDec(i,NonZero),NonZero);
            OffDec(i,~NonZero) = SBXhalf(Parent1(i,~NonZero),Parent2(i,~NonZero));
            OffDec(i,~NonZero) = PM(Problem,OffDec(i,~NonZero),~NonZero);
        end                       
    else
        OffDec = ones(size(OffMask));
    end
end

function Offspring = SBXhalf(Parent1,Parent2)
    [proC,disC,~,~] = deal(1,20,1,20);
    [N,D] = size(Parent1);
    beta  = zeros(N,D);
    mu    = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;
end

function Offspring = PM(Problem,Offspring,indexSite)
    [N,D] = size(Offspring);
    [proM,disM] = deal(1,20);
    % Polynomial mutation
    Lower = repmat(Problem.lower(indexSite),N,1);
    Upper = repmat(Problem.upper(indexSite),N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
 end