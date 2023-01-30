function FRR = CreditAssignment(SW,D)
% Credit assignment

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    K = 4;  % Number of operators
    Reward = zeros(1,K);
    for i = 1 : K
        Reward(i) = sum(SW(2,SW(1,:)==i));
    end
    [~,Rank] = sort(Reward,'descend');
    [~,Rank] = sort(Rank);
    Decay    = D.^Rank.*Reward;
    FRR      = Decay./sum(Decay);
end