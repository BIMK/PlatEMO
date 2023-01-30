function [Population,W,B] = updateWeight(Population,W,Z,T,EP,nus)
% Delete overcrowded subproblems and add new subproblems

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias
    
	%% Parameter setting
    [N,M] = size(Population.objs);    
    	
    %% Delete the overcrowded subproblems
    Dis = pdist2(Population.objs,Population.objs);
    Dis(logical(eye(length(Dis)))) = inf;
    Del = false(1,length(Population));
    while sum(Del) < min(nus,length(EP))
        Remain = find(~Del);
        subDis = sort(Dis(Remain,Remain),2);
        [~,worst] = min(prod(subDis(:,1:min(M,length(Remain))),2));
        Del(Remain(worst)) = true;
    end
    Population = Population(~Del);
    W = W(~Del,:);
	
    %% Add new subproblems
    % Determine the new solutions be added
    Combine  = [Population,EP];
    Selected = false(1,length(Combine));
    Selected(1:length(Population)) = true;
    Dis = pdist2(Combine.objs,Combine.objs);
    Dis(logical(eye(length(Dis)))) = inf;
    while sum(Selected) < min(N,length(Selected))
        subDis = sort(Dis(~Selected,Selected),2);
        [~,best] = max(prod(subDis(:,1:min(M,size(subDis,2))),2));
        Remain = find(~Selected);
        Selected(Remain(best)) = true;
    end
    % Add new subproblems
    newObjs = EP(Selected(length(Population)+1:end)).objs;
    temp    = 1./(newObjs-(repmat(Z,size(newObjs,1),1) - 0.001));
    temp(temp==inf) = 0.999999; %1e-6; % when (temp == Z) then 0
    W = [W;temp./repmat(sum(temp,2),1,size(temp,2))];
	
	% Redetect the neighbours of each solution
	B = pdist2(W,W);
	[~,B] = sort(B,2);
	B = B(:,1:T);
	
    % Add new solutions
    Population = Combine(Selected);
end