function [model,centers] = FCMmodel(Problem,Population,L1,L2)
% Fuzzy clustering-based method for modeling c_size* M models, where c_size
% is the number of clusters and M the number of objectives.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopDec      = Population.decs;
    csize       = 1 + ceil((length(Population)-L1)/L2);
    [centers,~] = fcm(PopDec,csize,[2 NaN 0.05 false]);
    dis         = pdist2(PopDec,centers);
    [~,index]   = sort(-dis);
    group       = index(1:L1,:);

    %% Build GP model of each objective for each cluster 
    model   = cell(csize,Problem.M);
    THETA   = 5.*ones(csize,Problem.M,Problem.D);
    for i   = 1 : csize
        temp = Population(group(:,i));
        PopDec = temp.decs;
        PopObj = temp.objs;
        for j = 1 : Problem.M
            dmodel = dacefit(PopDec,PopObj(:,j),...
                     'regpoly0','corrgauss',...
                     squeeze(THETA(i,j,:)),...
                     1e-5.*ones(1,Problem.D),...
                     100.*ones(1,Problem.D));
            model{i,j}   = dmodel;
            THETA(i,j,:) = dmodel.theta;
        end
    end
end