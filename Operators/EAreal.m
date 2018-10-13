function Offspring = EAreal(Global,Parent)
% <operator> <real>
% Simulated binary crossover and polynomial mutation
% proC ---  1 --- The probability of doing crossover
% disC --- 20 --- The distribution index of simulated binary crossover
% proM ---  1 --- The expectation of number of bits doing mutation 
% disM --- 20 --- The distribution index of polynomial mutation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [proC,disC,proM,disM] = Global.ParameterSet(1,20,1,20);
    Parent    = Parent([1:end,1:ceil(end/2)*2-end]);
    ParentDec = Parent.decs;
    [N,D]     = size(ParentDec);

    %% Simulated binary crossover
    Parent1Dec = ParentDec(1:N/2,:);
    Parent2Dec = ParentDec(N/2+1:end,:);
    beta = zeros(N/2,D);
    mu   = rand(N/2,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N/2,D);
    beta(rand(N/2,D)<0.5) = 1;
    beta(repmat(rand(N/2,1)>proC,1,D)) = 1;
    OffspringDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
                    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];

    %% Polynomial mutation
    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                         (1-(OffspringDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                         (1-(Upper(temp)-OffspringDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

    Offspring = INDIVIDUAL(OffspringDec);
end