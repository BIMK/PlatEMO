function Offspring = ARSBX(Problem,Parent,Parameter)
% Rotation based simulated binary crossover and polynomial mutation

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [B,Centroid,pc]       = deal(Parameter{:});
    [proC,disC,proM,disM] = deal(1,2,1,20);

    %% Simulated binary crossover
    ParentDec = Parent.decs;
    [N,D]     = size(ParentDec);  
    Norigin = round(N*pc/2)*2;
    Neig = round(N*(1-pc)/2)*2;
    ParentDecOrigin = ParentDec(1:Norigin,:);
    ParentDecEig = ParentDec(end-Neig+1:end,:);
    Parent1Dec = ParentDecOrigin(1:Norigin/2,:);
    Parent2Dec = ParentDecOrigin(Norigin/2+1:end,:);
    beta = zeros(Norigin/2,D);
    mu   = rand(Norigin/2,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],Norigin/2,D);
    beta(rand(Norigin/2,D)<0.5) = 1;
    beta(repmat(rand(Norigin/2,1)>proC,1,D)) = 1;
    OffspringDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
                    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];
    Flag = ones(Norigin,1);

    NormalParentDec = ParentDecEig - Centroid;
    [~,D] = size(B);
    ParentDecEig = NormalParentDec*B;
    Parent1Dec = ParentDecEig(1:Neig/2,:);
    Parent2Dec = ParentDecEig(Neig/2+1:end,:);
    beta = zeros(Neig/2,D);
    mu   = rand(Neig/2,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],Neig/2,D);
    beta(rand(Neig/2,D)<0.5) = 1;
    beta(repmat(rand(Neig/2,1)>proC,1,D)) = 1;
    OffspringRDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
                    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];
    OffspringRDec = OffspringRDec*B' + Centroid;
    Flag = [Flag;2*ones(Neig,1)];       
    OffspringDec = [OffspringDec;OffspringRDec];
    
    %% Polynomial mutation
    [~,D] = size(ParentDec);
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(nthroot(2.*mu(temp)+(1-2.*mu(temp)).*...
                         (1-(OffspringDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1),disM+1)-1);              
    temp = Site & mu>0.5; 
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(1-nthroot(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                        (1-(Upper(temp)-OffspringDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1),disM+1));
    Offspring = Problem.Evaluation(OffspringDec,Flag);
end