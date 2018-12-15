function [Population,Egamma] = EnvironmentalSelection(Population,Offspring,W,B,Z,gamma)
% The environmental selection of MOEA/D-MRDL

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Egamma = [];
    convergenceDirection = [];
    for i = randperm(length(Population))
        % Find the neighboring parents which have larger Tchebycheff scalar
        % function values than the offspring
        g_old = max(abs(Population(B(i,:)).objs-repmat(Z,size(B,2),1)).*W(B(i,:),:),[],2);
        g_new = max(abs(repmat(Offspring(i).obj-Z,size(B,2),1)).*W(B(i,:),:),[],2);
        PM    = B(i,g_old>g_new);
        % Replacement
        if ~isempty(PM)
            [~,nearest] = min(pdist2(Offspring(i).obj,Population(PM).objs));
            if isempty(convergenceDirection)
                % Record the convergence direction
                convergenceDirection = [convergenceDirection;Population(PM(nearest)).obj-Offspring(i).obj];
                % Replace the nearest parent with the offspring
                Population(PM(nearest)) = Offspring(i);
            else
                % Compute MRDL of each parent in PM to the offspring
                Sine1 = sqrt(1-(1-pdist2(Population(PM).objs,convergenceDirection,'cosine')).^2);
                Sine2 = sqrt(1-(1-pdist2(Offspring(i).obj,convergenceDirection,'cosine')).^2);
                RDL   = repmat(sqrt(sum(Population(PM).objs.^2,2)),1,size(Sine1,2)).*Sine1./repmat(norm(Offspring(i).obj).*Sine2,size(Sine1,1),1);
                MRDL  = max(RDL,[],2);
                if all(MRDL<=gamma)
                    % Record the convergence direction
                    convergenceDirection = [convergenceDirection;Population(PM(nearest)).obj-Offspring(i).obj];
                    % Replace the nearest parent with the offspring
                    Population(PM(nearest)) = Offspring(i);
                    % Record the MRDL for calculating the average MRDL
                    Egamma = [Egamma;MRDL(nearest)];
                end
            end
        end
    end
    % The average MRDL over the whole population
    Egamma = mean(Egamma);
end