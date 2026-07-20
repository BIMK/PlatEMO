function Archive = EnvironmentalSelection_SPEA2(Archive, N)
    ArchiveObj = Archive.objs;
    [ArchiveSize, ~] = size(Archive.objs);
    ideal_point      = min(ArchiveObj,[],1);
    worst_point      = max(ArchiveObj,[],1);   
    nf               = worst_point - ideal_point;
    nf(nf < 1e-100)  = 1e-100;
    normalized_obj   = ArchiveObj ./ nf;  
    Distance         = pdist2(normalized_obj, normalized_obj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,ArchiveSize);
    while sum(Del) < ArchiveSize - N
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
    Archive = Archive(~Del);
end