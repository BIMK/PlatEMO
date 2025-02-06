function [Population] = Final(Problem,Arch,N,W,score)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if length(Arch) <= N
        Population = Problem.Evaluation(Arch.decs.*Arch.masks);
    else
        RF = [];
        MF = [];
        while size(W,1) > 0
            [srch,W] = assign(Arch,W);
            [rf,mf,rch,wr] = Construction(srch,W,score);
            RF = [rf;RF];
            MF = [mf;MF];
            Arch(rch) = [];
            W(wr,:)   = [];
        end
        Population = Problem.Evaluation(RF.*MF);
    end
end

function [ar,w] = assign(arch,W)
    Popobj = arch.objs;
    L      = max(Popobj)-min(Popobj);
    objs   = (Popobj-min(Popobj))./L ;
    [~,n]  = size(arch);
    [N ,~] = size(W);
    for x = 1 : n
        obj = objs(x,:);
        for y = 1 : N
            s = sum(W(y,:).*obj,2);
            m = sqrt(sum(W(y,:).*W(y,:),2)*sum(obj.*obj,2));
            dang(1,y) = acos(s/m);
            [~,h]    = sort(dang);
        end
        arch(x).tno = h(1);
    end
    w  = W;
    ar = arch;
end