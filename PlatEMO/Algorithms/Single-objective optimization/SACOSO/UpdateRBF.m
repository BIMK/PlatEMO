function [tArchive,GbestRBF] = UpdateRBF(Problem,tArchive,SwarmRBF,GbestRBF)
% Update solutions in RBF-assisted swarm

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    [N,D] = size(SwarmRBF);
    D     = D - 1;
    [value,best] = min(SwarmRBF(:,1+D));
    
    new = Problem.Evaluation(SwarmRBF(best,1:D));
    SwarmRBF(best,D+1) = new.objs;
    if new.objs < GbestRBF(:,D+1)
        GbestRBF = SwarmRBF(best,:);
    end
    tArchive = [tArchive,new];
    
end