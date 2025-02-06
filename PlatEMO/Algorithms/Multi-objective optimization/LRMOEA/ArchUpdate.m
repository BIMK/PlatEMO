function [Arch] = ArchUpdate(Problem,arch,TemArch,thea,eta)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    TDec = arch.decs;
    Mask = arch.masks;
    % Perturb Solutions
    PopDec = TDec.*Mask;
    Population = Problem.Perturb(PopDec,1);

    Pop  = Population.objs;
    arch = Memory(arch,Pop);
    Mr   = arch.mrs;
    arch = arch(Mr<=thea);
    if isempty(arch)
        Arch = TemArch;
    else
        Tdec = (TemArch.decs).*(TemArch.masks);
        Decs = (arch.decs).*(arch.masks);
        L    = size(Tdec,1);
        n    = size(Decs,1);
        Tobj = TemArch.objs;
        a    = ones(L,1);
        for i = 1 : L
            for j = 1 : n
                if all(Tdec(i,:)==Decs(j,:))
                    a(i,1)  = 0;
                    arch(j) = Memory(arch(j),Tobj(i,:));
                    break;
                end
            end
        end
        t = TemArch(a==1);
        if ~isempty(t)
            maf   = max(t.objs);
            Tobjs = arch.objs;
            inx   = (sum(Tobjs,2)<(eta.*sum(maf)));
            arch  = arch(inx);
        end
        Arch = [arch,t];
    end
end

function value = Memory(arch,pop)
    for k = 1 : size(arch,2)
        arch(k).mobj = [arch(k).mobj;pop(k,:)];
        arch(k).Gn   = arch(k).Gn + 1;
        ob = arch(k).mobj;
        arch(k).mr = mean(mean(abs(ob-min(ob)))./mean(ob));
    end
    value = arch;
end