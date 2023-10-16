function Population = vssps(prob, sLower, sUpper)
% Randomly generate an initial population

% This function is written by Ian Meyer Kropp

    %  Nomenclature example
    %  N = 8 
    %  D = 14
    %  
    %          Cycle length of 14 
    %          |
    %  |-------|-----------------|
    %  1 1 1 1                        -
    %          1 1 1 1                |-- One full cycle 
    %                  1 1 1          |
    %                        1 1 1    -
    %     |-------|------|-----|---------------- Cycle count of 4
    %   |---|--|-|----------------------Cycle count of 4 
    %  1 1 
    %      1 1 
    %          1 
    %            1   
    %  |----|----|
    %       | 
    %       Cyle length of 6

    %% Result set up
    pop = prob.Initialization();
    varCount = size(prob.lower,2);
    mask = false(prob.N, varCount);

    %% Determine the positioning of each stripe per individual
    densityVector = 1 - linspace(sLower, sUpper, prob.N);

    widthVector = round(densityVector.*prob.D);
    
    % Put widths back into bound if rounding error occurred
    lb = floor((1- sLower)*prob.D);
    widthVector(widthVector > lb) = lb;

    cumulativeWidths = cumsum(widthVector);

    % if all sparsities are 100%, then skip processing, since everything
    % will be zeros 
    if sum(widthVector == 0) == prob.N
        processedIndvs = prob.N; 
    else
        processedIndvs = 0;
    end

    cycle_count = 0;
    cycles = zeros(prob.N, prob.D);
    while processedIndvs < prob.N
        % Figure out how many stripes will fit in this cycle 
        cycle_count = cycle_count + 1;
        spotsThatFitMask = cumulativeWidths <= prob.D & cumulativeWidths ~= 0;
        numThatFit = sum(spotsThatFitMask);
        largestFit = max(cumulativeWidths(spotsThatFitMask));

        cumulativeWidths = cumulativeWidths - largestFit; 
        cumulativeWidths(cumulativeWidths<0) = 0; 
        processedIndvs = processedIndvs + numThatFit;
        spotsThatFit = find(spotsThatFitMask);
        cycles(cycle_count,1:numThatFit) = spotsThatFit;   

    end

    %% Create density mask 

    % Mask out non-zero values cycle-by-cycle
    currentIndv = 1;
    for c = 1:cycle_count
        cycle = cycles(c, cycles(c, :) ~= 0);

        widths = widthVector(cycle);

        gapToFill = prob.D - sum(widths);
        gapSize = ceil((prob.D - sum(widths))/numel(widths));
        % For each individual in the cycle 
        position = 1;
        for i = 1:numel(widths)

            width = widths(i);

            % Determine if a gap is needed 
            gapWidth = 0;
            if gapToFill > 0
                gapWidth = gapSize;
                gapToFill = gapToFill - gapWidth;
            end

            % Determine the position of the stripe
            startPoint = position;

            if c == cycle_count
                endPoint = position+width-1; 
            else
                endPoint = position+width-1+gapWidth; 
            end

            % Prevent overflow from a gap calculation
            if endPoint > prob.D 
                endPoint = prob.D;
            end
    
            % Mask out stripe
            mask(currentIndv, startPoint:endPoint) = true; 

            % Go to the next individual 
            
            position = position+width+gapWidth;
            %position = position+width;

            currentIndv = currentIndv + 1;

        end

    end

    %% Mask off population according to stripe position 
    sparse_pop = pop.decs;
    sparse_pop(~mask) = 0;

    % Recalculate objective and constraints
    popDec = prob.CalDec(sparse_pop);
    popObj = prob.CalObj(sparse_pop);
    popCon = prob.CalCon(sparse_pop);

    Population = SOLUTION(popDec, popObj, popCon);
end