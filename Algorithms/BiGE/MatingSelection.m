function MatingPool = MatingSelection(BiObj)
% The mating selection of BiGE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(BiObj,1);
    
    %% Binary tournament selection
    Parents1   = randi(N,1,N);
    Parents2   = randi(N,1,N);
    Dominate   = any(BiObj(Parents1,:)<BiObj(Parents2,:),2) - any(BiObj(Parents1,:)>BiObj(Parents2,:),2);
    MatingPool = [Parents1(Dominate>=0),...
                  Parents2(Dominate<0)];
end