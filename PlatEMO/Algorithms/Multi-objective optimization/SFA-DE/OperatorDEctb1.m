function Offspring = OperatorDEctb1(Problem, Parent, best1, random1, random2, Parameter)
% The operator of differential evolution current-to-best/1

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yuma Horaguchi

    %% Parameter setting
    [CR, F] = deal(Parameter{:});
    [N, D]  = size(Parent);

    %% Differental evolution
    Site  = rand(N, D) <= CR;
    jrand = randi(D, [N, 1]);
    for i = 1 : N
        Site(i, jrand(i)) = 1;
    end
    Offspring       = Parent;
    Offspring(Site) = Parent(Site) + F * (best1(Site) - Parent(Site)) + F * (random1(Site) - random2(Site));
    Lower           = repmat(Problem.lower, N, 1);
    Upper           = repmat(Problem.upper, N, 1);
    Offspring       = min(max(Offspring, Lower), Upper);
end