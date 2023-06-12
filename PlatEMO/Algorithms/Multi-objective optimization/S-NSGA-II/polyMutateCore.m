function result = polyMutateCore(genome, lb, ub, eta)
% Ported from the polynomial mutation in Pymoo, circa early 2022

% This function is written by Ian Meyer Kropp

    delta1 = (genome - lb) / (ub - lb);    % Should be between 0 and 1
    delta2 = (ub - genome) / (ub - lb);    % Should be between 0 and 1

    exp = (eta + 1) ^ -1;
    
    ran = rand(size(genome));              
    deltaq = zeros(size(genome));
    
    leftMask = ran < 0.5;
    rightMask = ran >= 0.5;
    
    xy = 1 - delta1;
    val = 2.0 * ran + (1.0 - 2.0 * ran) .* (xy .^ (eta + 1.0));
    d = (val .^ exp) - 1.0;
    deltaq(leftMask) = d(leftMask);

    xy = 1.0 - delta2;
    val = 2.0 * (1.0 - ran) + 2.0 * (ran - 0.5) .* (xy .^ (eta + 1.0));
    d = 1.0 - (val .^ exp);
    deltaq(rightMask) = d(rightMask);
 
    muted_genome = genome + deltaq .* (ub - lb);

    muted_genome = min(max(muted_genome,lb),ub);
    
    result = muted_genome;
end