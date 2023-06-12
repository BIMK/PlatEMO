function newPop = sm2target(Pop, lb, ub, newSparsities)

% This function is written by Ian Meyer Kropp

    % Sparse Mutate to a target
    if numel(Pop) == 0 
        newPop = Pop;
        return;
    end

    [~,D] = size(Pop);

    sparsities = sum(Pop == 0, 2) / D;
    
    mutateMask = newSparsities ~= sparsities;

    indv2mut = find(mutateMask);

    % determine the non-zero increase/decrease for each indiv
    nz2add = round(D * (sparsities(mutateMask) - newSparsities(mutateMask)));
    
    % find where the non-zeros are
    [zIndvs,   zGenes] = find(Pop(indv2mut,:) == 0);
    % find where the zeros are 
    [nzIndvs, nzGenes] = find(Pop(indv2mut,:) ~= 0);
    
    newNzs = false(size(Pop));
    newZs = false(size(Pop));

    
    for i = 1:size(indv2mut,1)
        % gather relevant info
        indv_i = indv2mut(i);
        
        % Case where more non-zeros are needed
        if nz2add(i) > 0
            % find where there are zeros to flip
            zeroLocs = zGenes(zIndvs == i);
            
            % determine how many of them to flip
            numToFlip = nz2add(i);
            
            % Determine which of these posible zero positions to flip
            toFlip = zeroLocs(randperm(length(zeroLocs), numToFlip));

            % Record these positions 
            newNzs(indv_i, toFlip) = true;
            
        % Case where more zeros are needed
        else
            
            % find where there are non-zeros to flip
            nZeroLocs = nzGenes(nzIndvs == i);
            
            % determine how many of them to flip
            numToFlip = -nz2add(i);

            % Determine which of these posible non-zero positions to flip
            toFlip = nZeroLocs(randperm(length(nZeroLocs), numToFlip));

            % Record these positions 
            newZs(indv_i, toFlip) = true;

        end
    end

    %% Make the mutations
    % Find the min/max of the genome positions to mutate 
    [~, newNzsCols] = find(newNzs);

    newNzsLb = lb(newNzsCols);
    newNzsUb = ub(newNzsCols);

    Pop(newNzs) = newNzsLb + rand(1,sum(newNzs, 'all')) .* (newNzsUb - newNzsLb);
    Pop(newZs) = zeros(sum(newZs, 'all'), 1);

    newPop = Pop;
end