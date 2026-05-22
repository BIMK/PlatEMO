function p= Shape_Estimate(Population,N,z,znad)
% Estimate the  shape of PF

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [FrontNo,~] = NDSort(Population.objs,N);
    Population  = Population(FrontNo<=1);
    PopObj      = Population.objs ;
    [N,~]       = size(PopObj);
    PopObj      = (PopObj-repmat(z,size(PopObj,1),1))./repmat(znad-z ,size(PopObj,1),1);

    k  = 1.5;
    CP = 0.5 : 0.1 : 2;
    for i = 1 : length(CP)
        Gp   = (sum(PopObj.^CP(i),2)).^(1/CP(i));
        temp = sort(Gp);
        Q1   = temp(max(fix(N*0.25),1));
        Q3   = temp(max(fix(N*0.75),1));
        Max  = Q3+k*(Q3-Q1);
        Gp(Gp>Max) = []; % Gp is denoised using box plot
        Vp(i) = std(Gp./max(Gp));
    end
    [~,index] = min(Vp);
    p = CP(index);
end