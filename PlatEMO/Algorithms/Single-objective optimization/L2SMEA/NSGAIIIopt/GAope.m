function offDec = GAope(popDec,BU,BD)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    proC = 1;
    disC = 20;
    proM = 1;
    disM = 20;
    
    %% Select parents
    Parent1 = popDec(1:floor(end/2),:);
    Parent2 = popDec(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);
    
    %% Simulated binary crossover
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    offDec = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
              (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];

    %% Polynomial mutation
    Lower = repmat(BD,2*N,1);
    Upper = repmat(BU,2*N,1);
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    offDec       = min(max(offDec,Lower),Upper);
    offDec(temp) = offDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                              (1-(offDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    offDec(temp) = offDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                              (1-(Upper(temp)-offDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end