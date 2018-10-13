function [c,G] = PAES(Global,c,G,H,l_fails,l_opt,div)
% PAES local search

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    fails = 0;
    moves = 0;
    while fails < l_fails && moves < l_opt
        m = Global.Variation(c,1,@FEP);
        if all(c.obj<=m.obj)
            fails = fails + 1;
        else
            [H,dominated,~,mCrowd,cCrowd] = UpdateArchive(H,m,c,Global.N,div);
            if all(c.obj>=m.obj)
                c = m;
                fails = 0;
            elseif ~dominated && mCrowd < cCrowd
                c = m;
            end
        end
        G     = UpdateArchive(G,m,INDIVIDUAL(),Global.N,div);
        moves = moves + 1;
    end
end