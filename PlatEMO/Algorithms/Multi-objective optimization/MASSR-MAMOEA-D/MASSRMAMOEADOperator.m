function Offspring = MASSRMAMOEADOperator(Problem,current,parents,operator,F)
%MASSRMAMOEADOperator - Generate one MOEA/D differential offspring.
%
%   Choose one of four differential operators with fixed K=0.5 where
%   applicable. Repair bounds and apply polynomial mutation before
%   evaluating the offspring.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Differential mutation with the selected operator
    switch operator
        case 1
            mutant = current + F*(parents(1,:)-parents(2,:));
        case 2
            mutant = current + F*(parents(1,:)-parents(2,:)) + F*(parents(3,:)-parents(4,:));
        case 3
            mutant = current + 0.5*(current-parents(1,:)) + ...
                F*(parents(2,:)-parents(3,:)) + F*(parents(4,:)-parents(5,:));
        case 4
            mutant = current + 0.5*(current-parents(1,:)) + F*(parents(2,:)-parents(3,:));
    end

    % Bound repair and polynomial mutation
    lower  = reshape(Problem.lower,1,[]);
    upper  = reshape(Problem.upper,1,[]);
    mutant = min(max(mutant,lower),upper);
    span   = upper - lower;
    if isscalar(span)
        lower = repmat(lower,1,numel(current));
        upper = repmat(upper,1,numel(current));
        span  = upper - lower;
    end
    site          = rand(size(current)) < 1/numel(current) & span > 0;
    mu            = rand(size(current));
    left          = site & mu <= 0.5;
    delta         = (mutant(left)-lower(left))./span(left);
    mutant(left)  = mutant(left) + span(left).*((2*mu(left)+ ...
        (1-2*mu(left)).*(1-delta).^21).^(1/21)-1);
    right         = site & mu > 0.5;
    delta         = (upper(right)-mutant(right))./span(right);
    mutant(right) = mutant(right) + span(right).*(1-(2*(1-mu(right))+ ...
        2*(mu(right)-0.5).*(1-delta).^21).^(1/21));
    Offspring     = Problem.Evaluation(min(max(mutant,lower),upper));
end