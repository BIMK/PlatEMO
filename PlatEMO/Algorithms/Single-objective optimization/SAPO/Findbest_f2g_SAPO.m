function [bx, bid] = Findbest_f2g_SAPO(m, pop, fit, con)
%% Get the best solution
% 1st: gm <= 0 and arg min f(x)
% 2nd: gm >  0 and arg min gm(x) 

    ind = find(con(:, m) <= 0);
    if ~isempty(ind)    % 1st: gm <= 0 and arg min f(x)
        [~, bid] = min(fit(ind, 1));
        bx = pop(ind(bid), :);
    else                % 2nd: gm >  0 and arg min gm(x) 
        [~, bid] = min(con(:, m));
        bx = pop(bid, :);
    end
end