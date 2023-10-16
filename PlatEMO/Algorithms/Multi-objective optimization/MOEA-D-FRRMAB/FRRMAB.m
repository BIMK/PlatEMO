function op = FRRMAB(FRR,SW,C)
% Bandit-based operator selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if any(FRR==0) || any(SW(1,:)==0)
        op = randi(length(FRR));
    else
        n      = hist(SW(1,:),1:length(FRR));
        [~,op] = max(FRR+C*sqrt(2*log(sum(n))./n));
    end
end