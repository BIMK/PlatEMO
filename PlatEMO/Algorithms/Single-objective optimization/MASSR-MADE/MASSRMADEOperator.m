function Offspring = MASSRMADEOperator(Problem,Parents,operator,F,CR)
%MASSRMADEOperator - Generate one differential-evolution offspring.
%
%   operator selects rand/1, rand/2, current-to-best/1, or
%   current-to-rand/1. Apply binomial crossover with F and CR, then
%   polynomial mutation with probability 1/D and distribution index 20.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    x        = Parents(1).dec;
    P        = Parents(2:end).decs;
    [~,best] = min(FitnessSingle(Parents(2:end)));
    r        = randperm(size(P,1),5);

    % Differential mutation
    switch operator
        case 1
            v = P(r(1),:) + F*(P(r(2),:)-P(r(3),:));
        case 2
            v = P(r(1),:) + F*(P(r(2),:)-P(r(3),:)) + F*(P(r(4),:)-P(r(5),:));
        case 3
            v = x + F*(P(best,:)-x) + F*(P(r(1),:)-P(r(2),:));
        case 4
            v = x + rand*(P(r(1),:)-x) + F*(P(r(2),:)-P(r(3),:));
    end

    % Binomial crossover with at least one changed coordinate
    site                  = rand(size(x)) < CR;
    site(randi(numel(x))) = true;
    x(site)               = v(site);
    x                     = min(max(x,Problem.lower),Problem.upper);

    % Polynomial mutation and bound repair
    span      = Problem.upper - Problem.lower;
    site      = rand(size(x)) < 1/Problem.D & span > 0;
    u         = rand(size(x));
    left      = site & u <= 0.5;
    d         = (x(left)-Problem.lower(left))./span(left);
    x(left)   = x(left) + span(left).*((2*u(left)+(1-2*u(left)).*(1-d).^21).^(1/21)-1);
    right     = site & u > 0.5;
    d         = (Problem.upper(right)-x(right))./span(right);
    x(right)  = x(right) + span(right).*(1-(2*(1-u(right))+2*(u(right)-0.5).*(1-d).^21).^(1/21));
    Offspring = Problem.Evaluation(min(max(x,Problem.lower),Problem.upper));
end