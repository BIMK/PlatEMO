function Offspring = MOEADM2M_operator(Global,Parent)
% <operator> <real>
% Crossover and mutation used in MOEA/D-M2M
% proM --- 1 --- The expectation of number of bits doing mutation 

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    proM = Global.ParameterSet(1);
    Parent    = Parent([1:end,1:ceil(end/2)*2-end]);
    ParentDec = Parent.decs;
    [N,D]     = size(ParentDec);

    %% Crossover
    rc = (2*rand(N/2,1)-1).*(1-rand(N/2,1).^(-(1-Global.gen/Global.maxgen).^0.7));
    Parent1Dec   = ParentDec(1:N/2,:);
    Parent2Dec   = ParentDec(N/2+1:end,:);
    OffspringDec = Parent1Dec + repmat(rc,1,D).*(Parent1Dec-Parent2Dec);
    
    %% Mutation
    rm    = 0.25*(2*rand(N/2,D)-1).*(1-rand(N/2,D).^(-(1-Global.gen/Global.maxgen).^0.7));
    Site  = rand(N/2,D) < proM/D;
    Lower = repmat(Global.lower,N/2,1);
    Upper = repmat(Global.upper,N/2,1);
    OffspringDec(Site) = OffspringDec(Site) + rm(Site).*(Upper(Site)-Lower(Site));
                     
	%% Set the infeasible decision variables to feasible values
    temp1 = OffspringDec < Lower;
    temp2 = OffspringDec > Upper;
    rnd   = rand(N/2,D);
    OffspringDec(temp1) = Lower(temp1) + 0.5*rnd(temp1).*(Parent1Dec(temp1)-Lower(temp1));
    OffspringDec(temp2) = Upper(temp2) - 0.5*rnd(temp2).*(Upper(temp2)-Parent1Dec(temp2));

    Offspring = INDIVIDUAL(OffspringDec);
end