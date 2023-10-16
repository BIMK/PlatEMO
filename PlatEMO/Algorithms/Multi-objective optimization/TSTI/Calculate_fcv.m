function fcv = Calculate_fcv(Population)
% calculate normalized  constraints violation(CV) measuring feasibility

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    CV_Original = Population.cons;
    CV_Original(CV_Original<=0) = 0;
    CV = CV_Original./max(CV_Original);
    CV(:,isnan(CV(1,:))) = 0;
    fcv = sum(max(0,CV),2)./size(CV_Original,2);
end