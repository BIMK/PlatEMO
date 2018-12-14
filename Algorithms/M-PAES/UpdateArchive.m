function [G,dominated,GCrowd,oCrowd,pCrowd] = UpdateArchive(G,offspring,parents,N,div)
% Update the archive by the offspring, and return whether the offspring is
% dominated by the archive, the crowding degrees of solutions in the
% original archive, the crowding degree of the offspring, and the crowding
% degrees of the parents

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    oobj = repmat(offspring.obj,length(G),1);
    Gobj = G.objs;
    domi = any(oobj<=Gobj,2) - any(oobj>=Gobj,2);
    dominated              = any(domi==-1);
    [GCrowd,oCrowd,pCrowd] = GridDensity(div,Gobj,offspring.obj,parents.objs);
    if any(domi==1)
        G = [G(domi~=1),offspring];
    elseif ~dominated
        if length(G) < N
            G = [G,offspring];
        elseif any(oCrowd<GCrowd)
            worst    = find(GCrowd==max(GCrowd));
            worst    = worst(randi(length(worst)));
            G(worst) = offspring;
        end
    end
end