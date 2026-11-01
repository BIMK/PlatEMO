classdef CMOEAACRT< ALGORITHM
% <2026> <multi> <real> <constrained>
% Constrained multi-objective evolutionary algorithm with comprehensive ranking technique

%------------------------------- Reference --------------------------------
% Z. Zeng, Z. Wang, T. Wang, X. Z. Gao, S. Yang, and X. Yan. A constrained
% multi-objective evolutionary algorithm assisted by comprehensive ranking
% technique. Journal of Computational Design and Engineering, 2026, 13:
% 114–134
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)
    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            MP_pool    = [0.1;0.3;0.5;0.7;0.9];
            MP_index   = 1;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                t                   = Problem.FE/Problem.N;
                MaxGen              = Problem.maxFE/Problem.N;
                a                   = 0.5*(1-cos((1-t/MaxGen)*pi));
                [Offspring,M_recod] = Hybrid_DE(Population,Problem,MP_pool);
                [Population,~,~,result_sort] = Selection_ACRT([Population,Offspring],Problem.N,a);  % Environment selection method based on comprehensive ranking
                diff_sort    = result_sort(1:Problem.N)'-result_sort(Problem.N+1:2*Problem.N)';     % Improved individuals
                better_index = find(diff_sort>=1);  % obtain the inproved individuals
                M_recod_1    = M_recod;
                M_recod      = M_recod(better_index);
                M_uni        = unique(M_recod);
                if length(M_uni) >= 1
                    sel1      = length(find(M_recod==1))/length(find(M_recod_1==1));
                    sel2      = length(find(M_recod==2))/length(find(M_recod_1==2));
                    pro_sel_1 = sel1/(sel1+sel2);
                    pro_sel_1 = 0.4+0.5*pro_sel_1;  % Scaling the range of selection probabilities
                else
                    pro_sel_1 = 0.5;
                end
                MP_pool(MP_index) = pro_sel_1;
                MP_index          = MP_index+1;
                if MP_index > 5
                    MP_index = 1;
                end
            end
        end
    end
end