function [bx, bid] = Findbest_g2f_SAPO(pop, fit, cv)
%% Get the best solution
% Same to the Feasibility Rule

    ind = find(cv == 0);
    if ~isempty(ind)
        [~, bid] = min(fit(ind));
        bx = pop(ind(bid),:);
    else
        [~, bid] = min(cv);
        bx = pop(bid,:);
    end
end