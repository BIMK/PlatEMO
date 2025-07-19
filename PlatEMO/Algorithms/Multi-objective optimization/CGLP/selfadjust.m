function tipe = selfadjust(PopX,pop1,pop2,tipe)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    mkl1  = [];
    mkl2  = [];
    ttyh1 = pdist2(pop1,PopX');
    for ii = 1 : size(ttyh1,1)
        [mkl1(ii),~] = min(ttyh1(ii,:));
    end
    ggg1 = find(isnan(mkl1));
    if ~isempty(ggg1)
        mkl1(ggg1) = [];
        g1 = size(ttyh1,1)-size(ggg1,2);
    else
        g1 = size(ttyh1,1);
    end
    junzh1 = sum(mkl1)/g1;
    
    ttyh2 = pdist2(pop2,PopX');
    for ii = 1 : size(ttyh2,1)
        [mkl2(ii),~] = min(ttyh2(ii,:));
    end
    ggg2 = find(isnan(mkl2));
    if ~isempty(ggg2)
        mkl2(ggg2) = [];
        g2 = size(ttyh2,1)-size(ggg2,2);
    else
        g2 = size(ttyh2,1);
    end
    junzh2 = sum(mkl2)/g2;
    
    if junzh1 < junzh2
        if g2 > 3
            tipe = tipe+1;
        end
    else
        if g1>3
            tipe = tipe-1;
        end
    end
end