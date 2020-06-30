function MatingPool = MatingSelection(PopObj,Zmin)
% The mating selection of MaOEA-DDFC

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N  = size(PopObj,1);
    FC = CalFC(PopObj,Zmin);
    Parents1   = randi(N,1,N);
    Parents2   = randi(N,1,N);
    Dominate   = any(PopObj(Parents1,:)<PopObj(Parents2,:),2) - any(PopObj(Parents1,:)>PopObj(Parents2,:),2);
    MatingPool = [Parents1(Dominate==1),...
                  Parents2(Dominate==-1),...
                  Parents1(Dominate==0 & FC(Parents1)<=FC(Parents2)),...
                  Parents2(Dominate==0 & FC(Parents1)>FC(Parents2))];
end