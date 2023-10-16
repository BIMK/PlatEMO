function [Population,CrowdDis] = ArchiveUpdate(Population,N,eps,st)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenhua Li

    n = length(Population);

    if eps~=1 && st<0.5
       eps = 2*(1-eps)/(2*st+1)+2*eps-1;
    end

    %% Select global Pareto front
    [FrontNo,MaxFNo] = NDSort(Population.objs,n);
    next = FrontNo==1;
    first_pf = Population(next);
    new_pop = first_pf;
    remain_pop = Population(~next);

    V=0.2*prod(max(Population.decs)-min(Population.decs)).^(1./size(Population.decs,2));

    while ~isempty(remain_pop)
        %% Delete close solutions
        dist = min(pdist2(new_pop.decs,remain_pop.decs));
        index = dist<V;
        remain_pop(index) = [];
        if isempty(remain_pop)
            break;
        end

        %% Select remaining solutions
        [FrontNo,MaxFNo] = NDSort(remain_pop.objs,length(remain_pop));
        pick_pop = remain_pop(FrontNo==1);
        [nF,~] = NDSort([pick_pop.objs .* (1-eps); first_pf.objs],length(pick_pop)+length(first_pf));
        nF = nF(1:length(pick_pop));

        maxnF = max(nF);

        if maxnF>1
            new_pop = [new_pop pick_pop(nF==1)];
            remain_pop = remain_pop(FrontNo~=1);
            break;
        else
            new_pop = [new_pop pick_pop];
            remain_pop = remain_pop(FrontNo~=1);
        end
    end
    Population = new_pop;

    %% Balance the number of solutions in each Pareto front
    if length(Population) > N
        awd_index = [];
        [FrontNo,MaxFNo] = NDSort(Population.objs,length(Population));
        new_pop=[];
        n_sub_pop = ceil(N/MaxFNo);
        sel_pop = [];
        tmp_pop = [];
        for i=1:MaxFNo
            pop = Population(FrontNo==i);
            if length(pop)<n_sub_pop
                sel_pop = [sel_pop pop];
                awd_index = [awd_index (n_sub_pop-length(pop)).*ones(1,length(pop))];
            else
                tmp_pop = [tmp_pop pop];
            end
        end
        while length(tmp_pop) > N - length(sel_pop)
            dist = pdist2(tmp_pop.decs,tmp_pop.decs);
            dist = sort(dist);
            dist = sum(dist(1:3,:),1);
            [~,ind] = min(dist);
            tmp_pop(ind)=[];
        end
        awd_index = [awd_index zeros(1,length(tmp_pop))] + 1;
        Population = [sel_pop tmp_pop];
        CrowdDis = Crowding(Population.decs);
        CrowdDis = CrowdDis.* awd_index';
    else
        CrowdDis = Crowding(Population.decs);
    end
end