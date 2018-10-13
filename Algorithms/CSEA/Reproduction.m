function OffspringDec = Reproduction(Global,Ref,Next)
% Offspring repoduction by SBX and PM in CSEA

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He
    
    [proC,disC,proM,disM] = Global.ParameterSet(1,15,1,5);
    [~,D] = size(Next);

    %% Simulated binary crossover
	Decs = [Next;Ref.decs];
    Parent1Dec = Decs(1:floor(size(Decs,1)/2),:);
    Parent2Dec = Decs(floor(size(Decs,1)/2)+1:2*floor(size(Decs,1)/2),:);
    N   = floor(size(Decs,1)/2)*2;
    beta = zeros(N/2,D);
    miu  = rand(N/2,D);
    beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(disC+1));
    beta(miu>0.5)  = (2-2*miu(miu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N/2,D);
    beta(rand(N/2,D)<0.5) = 1;
    beta(repmat(rand(N/2,1)>proC,1,D)) = 1;
    OffspringDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
                    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];

    %% Polynomial mutation
    Site  = rand(N,D) < proM/D;
    miu   = rand(N,D);
    temp  = Site & miu<=0.5;
    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*((2.*miu(temp)+(1-2.*miu(temp)).*...
                         (1-(OffspringDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & miu>0.5; 
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-miu(temp))+2.*(miu(temp)-0.5).*...
                         (1-(Upper(temp)-OffspringDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end