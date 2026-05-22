function [PopDec,PopObj,MSE] = NSGAIII_opt(Achieve,Model,W,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

    PopDec = Achieve.decs;
    PopObj = Achieve.objs;
    Zmin   = min(PopObj,[],1);
    if size(PopDec,1) >= Problem.N
        Next   = nsga3EnvironmentalSelection(PopDec,PopObj,Problem.N,W,Zmin);
        PopDec = PopDec(Next,:);
        PopObj = PopObj(Next,:);
    end
    g    = 1;
    gmax = 20;
    beta = 0;
    while g <= gmax
        OffDec = generateOffing(PopDec,Problem,Achieve);
        PopDec = cat(1,PopDec,OffDec);
        [PopObj,~,MSE,~] = GP_estimate(PopDec,Model,Problem.M,beta);
        Zmin  = min([Zmin; PopObj],[],1);
        Choose = nsga3EnvironmentalSelection(PopDec,PopObj,Problem.N,W,Zmin);
        PopDec = PopDec(Choose,:);
        PopObj = PopObj(Choose,:);
        MSE    = MSE(Choose,:);
        g      = g + 1;
    end
end