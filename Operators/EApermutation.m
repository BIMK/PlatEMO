function Offspring = EApermutation(Global,Parent)
% <operator> <permutation>
% Order crossover and slight mutation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Parent    = Parent([1:end,1:ceil(end/2)*2-end]);
    ParentDec = Parent.decs;
    [N,D]     = size(ParentDec);
    
    %% Order crossover
    OffspringDec = ParentDec;
    k = randi(D,1,N);
    for i = 1 : N/2
        OffspringDec(i,k(i)+1:end)     = setdiff(ParentDec(i+N/2,:),ParentDec(i,1:k(i)),'stable');
        OffspringDec(i+N/2,k(i)+1:end) = setdiff(ParentDec(i,:),ParentDec(i+N/2,1:k(i)),'stable');
    end
    
    %% Slight mutation
    k = randi(D,1,N);
    s = randi(D,1,N);
    for i = 1 : N
        if s(i) < k(i)
            OffspringDec(i,:) = OffspringDec(i,[1:s(i)-1,k(i),s(i):k(i)-1,k(i)+1:end]);
        elseif s(i) > k(i)
            OffspringDec(i,:) = OffspringDec(i,[1:k(i)-1,k(i)+1:s(i)-1,k(i),s(i):end]);
        end
    end
    
    Offspring = INDIVIDUAL(OffspringDec);
end