function [FEA,INFEA] = Variable_classification(Problem,Population,C)
% Detect the kind of each decision variable

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    pop_dec = Population.decs;
    D       = size(pop_dec,2);
    PN      = 4;
    SN      = 5;

    theta = 0.00001;
    for d = 1 : D
        Per = (Problem.lower(d):(Problem.upper(d)-Problem.lower(d))/(PN-1):Problem.upper(d))';
        for j = 1 : SN
            pop_dec_SN      = repmat(C(j,:),PN,1);
            pop_dec_SN(:,d) = Per;
            pop_fenxi       = Problem.Evaluation(pop_dec_SN);
            pop_fenxi_con   = pop_fenxi.cons;
            pop_fenxi_con(pop_fenxi_con < 0) = 0;
            pop_fenxi_con   = sum(pop_fenxi_con,2);
            pop_dec_con     = pop_fenxi_con;
            VarCon(j,d)     = std(pop_dec_con);
        end
    end
    meanCon = mean(VarCon);
    FEA     = find(meanCon>theta);
    INFEA   = find(meanCon<=theta);
end