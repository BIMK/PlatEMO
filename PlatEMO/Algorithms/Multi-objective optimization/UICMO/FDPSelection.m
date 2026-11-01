function Population = FDPSelection(OldPop,NewPop,W,Z,ReferenceSet,current)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    N                = length(OldPop);
    AngleNew         = acos(1-pdist2(NewPop.objs,W,'cosine'));
    [~,AssociateNew] = min(AngleNew,[],2);
    AngleOld         = acos(1-pdist2(OldPop.objs,W,'cosine'));
    [~,AssociateOld] = min(AngleOld,[],2);
    CVOld            = sum(max(0,OldPop.cons),2);
    CVNew            = sum(max(0,NewPop.cons),2);
    Weights          = W(current,:);
    FitnessOld       = max(abs(OldPop.objs-repmat(Z,N,1)).*Weights,[],2);
    FitnessNew       = max(abs(NewPop.objs-repmat(Z,1,1)).*Weights,[],2);
    Population       = OldPop;
    replacements     = 0;
    for i = 1 : N
        if replacements == 2
            break;
        end
        better = false;
        if AssociateOld(i) == current(i) && AssociateNew == current(i)
            RegionOld = ReferenceRegion(OldPop(i),ReferenceSet);
            RegionNew = ReferenceRegion(NewPop,ReferenceSet);
            if RegionNew == 0 && RegionOld ~= 0
                better = true;
            elseif (RegionNew == 0 && RegionOld == 0) || ...
                    (RegionNew == 1 && RegionOld ~= 0) || ...
                    (RegionNew == -1 && RegionOld == 1)
                better = CVOld(i) > CVNew || ...
                    (CVOld(i) == CVNew && FitnessOld(i) > FitnessNew(i));
            elseif RegionNew == -1 && RegionOld == -1
                better = FitnessOld(i) > FitnessNew(i);
            end
        elseif AngleOld(i,current(i)) > AngleNew(current(i)) && FitnessOld(i) > FitnessNew(i)
            better = true;
        end
        if better
            Population(i) = NewPop;
            replacements  = replacements + 1;
        end
    end
end

function K = ReferenceRegion(Pop,ReferenceSet)
    N = size(ReferenceSet,2);
    K = 0;
    if ~isempty(ReferenceSet)
        for i = 1 : N
            k = any(Pop.objs<ReferenceSet(i).objs) - any(Pop.objs>ReferenceSet(i).objs);
            if k == 1
                K = 1;
                break
            elseif k == -1
                K = -1;
                break
            end
        end
    end
end