function Offspring = Operator(Parent1,Parent2,Parent3,Parameter)
% Differental evolution with variable-wise mutation restriction

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

    %% Parameter setting
    if nargin > 3
        [CR,proM,disM] = deal(Parameter{:});
    else
        [CR,proM,disM] = deal(0.15,1,20);
    end
    Parent1 = Parent1.decs;
    Parent2 = Parent2.decs;
    Parent3 = Parent3.decs;
    [N,D]   = size(Parent1);
    Global  = GLOBAL.GetObj();

    %% Calculate the difference vectors
    V = Parent2 - Parent3;
    
    %% Do polynomial mutation on the difference vectors
    Lower   = repmat(Global.lower,N,1);
    Upper   = repmat(Global.upper,N,1);
    Site    = rand(N,D) < proM/D;
    mu      = rand(N,D);
    temp    = Site & mu<=0.5;
    V(temp) = V(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)).^(1/(disM+1))-1);
    temp    = Site & mu>0.5; 
    V(temp) = V(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))).^(1/(disM+1)));

    %% Restrict the difference vectors
    V = min(max(V,(Lower-Upper)/2),(Upper-Lower)/2);

    %% Generate offsprings
    Site         = rand(N,D) < CR;
    OffDec       = Parent1;
    OffDec(Site) = OffDec(Site) + V(Site);
    Offspring    = INDIVIDUAL(OffDec);
end