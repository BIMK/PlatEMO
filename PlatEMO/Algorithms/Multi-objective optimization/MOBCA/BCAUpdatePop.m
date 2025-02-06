function new_armies = BCAUpdatePop(Archive,Population,nArmies,div)
% Update the archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Find the non-dominated solutions
    newpop=[Archive,Population];
    ndpop=newpop(NDSort(newpop.objs,1)==1);
    %% Grid-based retention
    if length(ndpop) > nArmies
        Del = Delete(ndpop.objs,length(ndpop)-nArmies,div);
        ndpop(Del) = [];
    end
    if size(ndpop,2)>=nArmies
        new_armies=ndpop(1:nArmies);
    else
%         a=[ndpop,Population];
%         [~,loc]=unique(a.decs,"rows");
%         new_armies=a(loc(1:nArmies));
        replace_flag=rand(1,nArmies)<0.1;
        replace_size=size(find(replace_flag==1),2);
        min_size=min(size(ndpop,2),replace_size);

        replace_index=find(replace_flag==1);
        index_index=randperm(min_size,min_size);
        Population(replace_index(index_index))=ndpop(1:min_size);
        new_armies=Population;
    end
end

function Del = Delete(PopObj,K,div)   
    N = size(PopObj,1);

    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);
    fmin = min(PopObj,[],1);
    d    = (fmax-fmin)/div;
    GLoc = floor((PopObj-repmat(fmin,N,1))./repmat(d,N,1));
    GLoc(GLoc>=div)   = div - 1;
    GLoc(isnan(GLoc)) = 0;

    %% Calculate the crowding degree of each grid
    [~,~,Site] = unique(GLoc,'rows');
    CrowdG     = hist(Site,1:max(Site));

    %% Delete K solutions
    Del = false(1,N);
    while sum(Del) < K
        % Select the most crowded grid
        maxGrid = find(CrowdG==max(CrowdG));
        Temp    = randi(length(maxGrid));
        Grid    = maxGrid(Temp);
        % And delete one solution randomly from the grid
        InGrid  = find(Site==Grid);
        Temp    = randi([1,length(InGrid)]);
        p       = InGrid(Temp);
        Del(p)  = true;
        Site(p) = NaN;
        CrowdG(Grid) = CrowdG(Grid) - 1;
    end
end