function Offspring = EAbinary(Global,Parent)
% <operator> <binary>
% One point crossover and bitwise mutation
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
    
    %% One point crossover
    Parent1Dec = ParentDec(1:N/2,:);
    Parent2Dec = ParentDec(N/2+1:end,:);
    k = repmat(1:D,N/2,1) > repmat(randi(D,N/2,1),1,D);
    Offspring1Dec    = Parent1Dec;
    Offspring2Dec    = Parent2Dec;
    Offspring1Dec(k) = Parent2Dec(k);
    Offspring2Dec(k) = Parent1Dec(k);
    OffspringDec     = [Offspring1Dec;Offspring2Dec];
    
    %% Bitwise mutation
    Site = rand(N,D) < proM/D;
    OffspringDec(Site) = ~OffspringDec(Site);
    
    Offspring = INDIVIDUAL(OffspringDec);
end