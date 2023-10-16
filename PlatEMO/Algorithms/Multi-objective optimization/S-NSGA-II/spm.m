function newPop = spm(Pop, lb, ub, Parameter)

% This function is written by Ian Meyer Kropp

    % Each row is a different population member
    % Each column is a different genome
    
    if nargin > 3
        [probMut,distrMut, probSMut, distrSMut] = deal(Parameter{:});
    else
        [probMut,distrMut, probSMut, distrSMut] = deal(1,20,1,20);
    end

    [N,D] = size(Pop);
    
    % Determine where the zeros are
    nonZeroMask = Pop ~= 0;
   
    %% Value mutations 
    ran = rand(size(Pop(nonZeroMask)));

    toMutateNZ = ran < (probMut/D);
    
    toMutate = false(size(Pop));
    
    toMutate(nonZeroMask) = toMutateNZ;
    
    [~, genomesToMutate] = find(toMutate);

    lb_pm = lb(genomesToMutate)';
    ub_pm = ub(genomesToMutate)';

    % mutate values
    Pop(toMutate) = polyMutateCore(  Pop(toMutate), ...
                                    lb_pm, ub_pm, distrMut);
        
    %% Sparsity mutations 

    % Determine which population members to mutate sparsity 
    ran = rand(N, 1);

    mutateMask = ran < probSMut/D;

    % Figure out the individual sparsities of each individual 
    sparsities = sum(Pop == 0, 2) / D;
    
    lb_sp = zeros(sum(mutateMask), 1);
    ub_sp = ones(sum(mutateMask), 1);

    newSparsities = sparsities;

    newSparsities(mutateMask) = polyMutateCore(sparsities(mutateMask), lb_sp, ub_sp, distrSMut);

    newSparsities = min(max(newSparsities,0),1);

    % check if there's anything to do
    if newSparsities == sparsities
        newPop = Pop;
    else
        newPop = sm2target(Pop, lb, ub, newSparsities);
    end
end