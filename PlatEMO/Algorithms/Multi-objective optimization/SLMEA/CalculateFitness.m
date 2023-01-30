function [Fitness,c,Z,O] = CalculateFitness(ArchiveMask,D)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    len     = length(ArchiveMask(:,1));
    ratio   = sum(ArchiveMask,1)./len;
    R       = abs(ratio-0.5);
    [a,~]   = min(R,[],2);
    EL      = find(R==a);
    ll      = length(EL);
    lo      = randperm(ll,1);
    vector  = ArchiveMask(:,EL(lo));
    vector  = repmat(vector,1,D);
    Hd      = sum(xor(vector,ArchiveMask),1);
    Ad      = sum(and(vector,ArchiveMask),1);
    Fitness = Hd./(Hd+Ad);
    c       = find(ratio>0&ratio<1);
    Z       = find(ratio==0);
    O       = find(ratio==1);
end