function Offspring = Operator(Parent1,Parent2)
% Crossover and mutation used in MOEA/D-M2M
% proM --- 1 --- The expectation of number of bits doing mutation 

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Parent1 = Parent1.decs;
    Parent2 = Parent2.decs;
    [N,D]   = size(Parent1);
    Global  = GLOBAL.GetObj();

    %% Crossover
    rc = (2*rand(N,1)-1).*(1-rand(N,1).^(-(1-Global.gen/Global.maxgen).^0.7));
    OffDec = Parent1 + repmat(rc,1,D).*(Parent1-Parent2);
    
    %% Mutation
    rm    = 0.25*(2*rand(N,D)-1).*(1-rand(N,D).^(-(1-Global.gen/Global.maxgen).^0.7));
    Site  = rand(N,D) < 1/D;
    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    OffDec(Site) = OffDec(Site) + rm(Site).*(Upper(Site)-Lower(Site));
                     
	%% Set the infeasible decision variables to feasible values
    temp1 = OffDec < Lower;
    temp2 = OffDec > Upper;
    rnd   = rand(N,D);
    OffDec(temp1) = Lower(temp1) + 0.5*rnd(temp1).*(Parent1(temp1)-Lower(temp1));
    OffDec(temp2) = Upper(temp2) - 0.5*rnd(temp2).*(Upper(temp2)-Parent1(temp2));
    Offspring     = INDIVIDUAL(OffDec);
end