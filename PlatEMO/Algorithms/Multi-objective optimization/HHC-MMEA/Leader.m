function leader = Leader(Mask,Front,GV,D)
% update the leader by GV

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    fs= zeros(1,D);
    sample = rand(1,D);
    fs(GV>sample) = 1;
    mask = Mask(Front==1,:);
    dis = pdist2(mask,fs,'hamming');
    index = find(dis==min(dis),1);
    leader = mask(index,:);
end