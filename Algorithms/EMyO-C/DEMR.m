function Offspring = DEMR(Global,Parent)
% <operator> <real>
% Differental evolution with variable-wise mutation restriction
% CR   --- 0.15 --- Parameter CR in differental evolution
% proM ---    1 --- The expectation of number of bits doing mutation 
% disM ---   20 --- The distribution index of polynomial mutation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

    [CR,proM,disM] = Global.ParameterSet(0.15,1,20);
    Parent    = Parent([1:end,1:ceil(end/3)*3-end]);
    ParentDec = Parent.decs;
    [N,D]     = size(ParentDec);

    %% Calculate the difference vectors
    Parent1Dec = ParentDec(1:N/3,:);
    Parent2Dec = ParentDec(N/3+1:N/3*2,:);
    Parent3Dec = ParentDec(N/3*2+1:end,:);
    V          = Parent2Dec - Parent3Dec;
    
    %% Do polynomial mutation on the difference vectors
    Lower   = repmat(Global.lower,N/3,1);
    Upper   = repmat(Global.upper,N/3,1);
    Site    = rand(N/3,D) < proM/D;
    mu      = rand(N/3,D);
    temp    = Site & mu<=0.5;
    V(temp) = V(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)).^(1/(disM+1))-1);
    temp    = Site & mu>0.5; 
    V(temp) = V(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))).^(1/(disM+1)));

    %% Restrict the difference vectors
    V = min(max(V,(Lower-Upper)/2),(Upper-Lower)/2);

    %% Generate offsprings
    Site               = rand(N/3,D) < CR;
    OffspringDec       = Parent1Dec;
    OffspringDec(Site) = OffspringDec(Site) + V(Site);
    
    Offspring = INDIVIDUAL(OffspringDec);
end