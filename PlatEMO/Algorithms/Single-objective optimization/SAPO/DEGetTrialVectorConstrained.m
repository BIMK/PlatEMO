function [u_pop] = DEGetTrialVectorConstrained(pop, fit, cv, F, CR, mut, xov, bid)
%% DE gets trial vectors

    [N,D] = size(pop);
    if isempty(bid)  % No input of bid
        % Feasibility Rule
        ind = find(cv == 0);
        if ~isempty(ind)
            [~, bid] = min(fit(ind));
            bx = pop(ind(bid),:);
        else
            [~, bid] = min(cv);
            bx = pop(bid,:);
        end
    else    % No need to select the best
        bx = pop(bid,:);
    end
    switch mut
        case 1  % rand/1
            indices = zeros(N, 3);
            for i = 1 : N
                indices(i, :) = datasample([1:i-1 i+1:N], 3, 'Replace', false);
            end
            v_pop = pop(indices(:, 1),:) + F.*(pop(indices(:, 2),:) - pop(indices(:, 3),:));
        case 2  % rand/2
            indices = zeros(N, 5);
            for i = 1 : N
                indices(i, :) = datasample([1:i-1 i+1:N], 5, 'Replace', false);
            end
            v_pop = pop(indices(:, 1),:) + F.*(pop(indices(:, 2),:) - pop(indices(:, 3),:)) + F.*(pop(indices(:, 4),:) - pop(indices(:, 5),:));
        case 3  % best/1
            indices = zeros(N, 2);
            for i = 1 : bid
                indices(i, :) = datasample([1:i-1 i+1:bid-1 bid+1:N], 2, 'Replace', false);
            end
            for i = bid+1 : N
                indices(i, :) = datasample([1:bid-1 bid+1:i-1 i+1:N], 2, 'Replace', false);
            end
            v_pop = repmat(bx, N, 1) + F.*(pop(indices(:, 1),:) - pop(indices(:, 2),:));
        case 4  % best/2
            indices = zeros(N, 4);
            for i = 1 : bid
                indices(i, :) = datasample([1:i-1 i+1:bid-1 bid+1:N], 4, 'Replace', false);
            end
            for i = bid+1 : N
                indices(i, :) = datasample([1:bid-1 bid+1:i-1 i+1:N], 4, 'Replace', false);
            end
            v_pop = repmat(bx, N, 1) + F.*(pop(indices(:, 1),:) - pop(indices(:, 2),:)) + F.*(pop(indices(:, 3),:) - pop(indices(:, 4),:));
        case 5  % current-to-rand/1
            indices = zeros(N, 3);
            for i = 1 : N
                indices(i, :) = datasample([1:i-1 i+1:N], 3, 'Replace', false);
            end
            v_pop = pop + F.*(pop(indices(:, 1),:) - pop) + F.*(pop(indices(:, 2),:) - pop(indices(:, 3),:));
        case 6  % current-to-best/1
            indices = zeros(N, 2);
            for i = 1 : bid
                indices(i, :) = datasample([1:i-1 i+1:bid-1 bid+1:N], 2, 'Replace', false);
            end
            for i = bid+1 : N
                indices(i, :) = datasample([1:bid-1 bid+1:i-1 i+1:N], 2, 'Replace', false);
            end
            v_pop = pop + F.*(repmat(bx, N, 1) - pop) + F.*(pop(indices(:, 1),:) - pop(indices(:, 2),:));
        case 7 % current-to-pbest/1
        case 8 % rand-to-best/1
            indices = zeros(N, 3);
            for i = 1 : bid
                indices(i, :) = datasample([1:i-1 i+1:bid-1 bid+1:N], 3, 'Replace', false);
            end
            for i = bid+1 : N
                indices(i, :) = datasample([1:bid-1 bid+1:i-1 i+1:N], 3, 'Replace', false);
            end
            v_pop = pop(indices(:, 1),:) + F.*(repmat(bx, N, 1) - pop(indices(:, 1),:)) + F.*(pop(indices(:, 2),:) - pop(indices(:, 3),:));
    end
    switch xov
        case 0  % no crossover
        case 1  % binomial
            mask  = rand(N, D) <= CR;
            randj = randi(D, N, 1);
            mask((1:N).' + (randj - 1) * N) = true;
            u_pop = v_pop .* mask + pop .* ~mask;
        case 2  % exponential
    end
end