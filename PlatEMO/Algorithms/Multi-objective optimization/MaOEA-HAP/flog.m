function k = flog(DA,N,p)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Popobj = DA.objs;
    M      = size(Popobj,2);
    
    %% Normalization
    Zmin     = min(Popobj,[],1);
    Zmax     = max(Popobj,[],1);
    Lp       = Shape_Estimate(DA,N,Zmin,Zmax);
    PopobjLP = DA.objs - Zmin;
    PopobjLP = PopobjLP+10^-6;
    tran_Obj = PopobjLP./repmat((sum(PopobjLP.^Lp,2)).^(1/Lp),1,M);
    for i = 1 : size(PopobjLP,1)
        DAPT(i) = norm(PopobjLP(i,:) - tran_Obj(i,:),p);
    end
    if all(DAPT<1)
        k = 1;
    else
        k = 0;
    end
end