function OffDec = generateOffing(PopDec,Problem,A1)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

    if size(PopDec,1) < Problem.N
        MatingPool = randi(size(PopDec,1),1,Problem.N);
        OffDec     = OperatorGA(Problem,PopDec(MatingPool,:));
    else
        OffDec = OperatorGA(Problem,PopDec);
    end
    pop_candi = [];
    NP        = size(OffDec,1);
    for ii = 1 : NP
        if min(sqrt(sum((OffDec(ii,:) - [A1.decs;pop_candi]).^2,2)))>1E-6
            pop_candi = cat(1,pop_candi,OffDec(ii,:));
        end
    end
    OffDec = pop_candi;
end