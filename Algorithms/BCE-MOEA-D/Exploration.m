function Offspring = Exploration(Global,PC,NPC,nND)
% Individual exploration in BCE

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PCObj  = PC.objs;
    NPCObj = NPC.objs;
    
    %% Normalization
    fmax   = max(PCObj,[],1);
    fmin   = min(PCObj,[],1);
    PCObj  = (PCObj-repmat(fmin,size(PCObj,1),1))./repmat(fmax-fmin,size(PCObj,1),1);
    NPCObj = (NPCObj-repmat(fmin,size(NPCObj,1),1))./repmat(fmax-fmin,size(NPCObj,1),1);

    %% Determine the size of the niche
    d  = pdist2(PCObj,PCObj);
    d(logical(eye(length(d)))) = inf;
    d  = sort(d,2);
    r0 = mean(d(:,min(3,size(d,2))));
    r  = nND/Global.N*r0;
    
    %% Detect the solutions in PC to be explored
    d = pdist2(PCObj,NPCObj);
    S = find(sum(d<=r,2)<=1);
    
    %% Generate new solutions
    if ~isempty(S)
        MatingPool = randi(length(PC),1,length(S));
        Offspring  = Global.Variation(PC([S',MatingPool]),length(S));
    else
        Offspring = [];
    end
end