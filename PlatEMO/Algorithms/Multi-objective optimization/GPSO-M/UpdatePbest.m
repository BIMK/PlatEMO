function Pbest = UpdatePbest(Pbest,Population)
% Update the local best position of each particle

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    temp = Pbest.objs - Population.objs;
    Dominate = any(temp<0,2) - any(temp>0,2);
    Pbest(Dominate==-1) = Population(Dominate==-1);
    temp = rand(length(Dominate),1);
    Pbest(Dominate==0 & temp<0.5) = Population(Dominate==0 & temp<0.5);
end