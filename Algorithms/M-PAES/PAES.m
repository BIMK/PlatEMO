function [c,G] = PAES(c,G,H,N,l_fails,l_opt,div)
% PAES local search

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    fails = 0;
    moves = 0;
    while fails < l_fails && moves < l_opt
        m = FEP(c);
        if all(c.obj<=m.obj)
            fails = fails + 1;
        else
            [H,dominated,~,mCrowd,cCrowd] = UpdateArchive(H,m,c,N,div);
            if all(c.obj>=m.obj)
                c = m;
                fails = 0;
            elseif ~dominated && mCrowd < cCrowd
                c = m;
            end
        end
        G     = UpdateArchive(G,m,INDIVIDUAL(),N,div);
        moves = moves + 1;
    end
end