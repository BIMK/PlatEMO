function [Z, Population, Archive, Sample, SVM] = initialize(Global)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

global incremental_sequence stable_threshold LARGE_ARCHIVE_FLAG SVM_KERNEL SVM_KERNEL_SCALE SVM_EXIST_FLAG reserved_evaluations reference_generation_mode crowding_pick_flag lb ub reference_generation_intial_flag learning_initial_flag normalization_str status_initial_flag disable_flag hp_initial_flag;
rng('shuffle');
SVM_EXIST_FLAG = false; SVM_KERNEL = 'gaussian'; SVM_KERNEL_SCALE = 0.056;
lb = zeros(1, Global.M); ub = ones(1, Global.M);
crowding_pick_flag = false;
learning_initial_flag = true;
reference_generation_intial_flag = true;
status_initial_flag = true;
hp_initial_flag = true;
reference_generation_mode = 'strict';
LARGE_ARCHIVE_FLAG = false;
reserved_evaluations = inf;
disable_flag = false;
SVM.a = [];
SVM.b = [];
SVM.g = [];
SVM.ind = [];
SVM.uind = [];
SVM.X_mer = [];
SVM.y_mer = [];
SVM.Rs = [];
SVM.Q = [];
SVM.scale = 0.056; SVM.type = 5; % 5 for Gaussian!
SVM.C = 10;
% if Global.M > 10
%     error('M > 10 not yet supported on PlatEMO_v2');
% end
problem = class(Global);
%% hyperparameters for MaF
% if ~regexp(problem, 'MaF')
%     error('problems other than MaF not yet supported on PlatEMO_v2');
% else
    if strcmp(problem, 'MaF1')
        if norm(stable_threshold) < 1e-2
            stable_threshold = [20, 20, 15];
        end
        if Global.M == 5
            disable_flag = false;
            crowding_pick_flag = false;
            incremental_sequence = [210,330,495,715,1001,1365,1820,2380,3060,3876,4845,8855];
        elseif Global.M == 10
            incremental_sequence = [230,275,440,16445];
            disable_flag = false;
            crowding_pick_flag = true;
        end
        normalization_str = '';
        reserved_evaluations = 150000;
    end
    if strcmp(problem, 'MaF2')
        if Global.M == 5
            incremental_sequence = [210;330;495;715;1001;1365;1820;2380];
            disable_flag = false;
            crowding_pick_flag = false;
        elseif Global.M == 10
            incremental_sequence = [230,275,440,715,725,770,935,1430,2002];
            disable_flag = false;
            crowding_pick_flag = true;
            normalization_str = '';
        end
    end
    if strcmp(problem, 'MaF3')
        if norm(stable_threshold) < 1e-2
            stable_threshold = [inf, inf, inf];
        end
        crowding_pick_flag = false;
        disable_flag = true;
        normalization_str = '';
    end
    if strcmp(problem, 'MaF4')
        if norm(stable_threshold) < 1e-2
            stable_threshold = [5, 10, 15];
        end
        if Global.M == 5
            disable_flag = false;
            LARGE_ARCHIVE_FLAG = true;
            crowding_pick_flag = false;
        elseif Global.M == 10
            disable_flag = false;
            crowding_pick_flag = false;
            reserved_evaluations = 150000;
        end
        normalization_str = '';
    end
    if strcmp(problem, 'MaF5')
        crowding_pick_flag = false;
        disable_flag = true;
        normalization_str = 'normalize';
    end
    if strcmp(problem, 'MaF6')
        if norm(stable_threshold) < 1e-2
            stable_threshold = [5, 10, 15];
        end
        if Global.M == 5
            disable_flag = false;
            crowding_pick_flag = true;
        elseif Global.M == 10
            disable_flag = false;
            crowding_pick_flag = false;
        end
        normalization_str = '';
        reference_generation_mode = 'normal';
    end
    if strcmp(problem, 'MaF7')
        if norm(stable_threshold) < 1e-2
            stable_threshold = [10, 15, 20];
        end
        if Global.M == 5
            disable_flag = false;
            crowding_pick_flag = false;
        elseif Global.M == 10
            disable_flag = false;
            crowding_pick_flag = false;
        end
        normalization_str = 'normalize';
        reserved_evaluations = 150000;
        reference_generation_mode = 'normal';
    end
    if strcmp(problem, 'MaF8')
        if norm(stable_threshold) < 1e-2
            stable_threshold = [5, 5, 5];
        end
        if Global.M == 5
            disable_flag = false;
            crowding_pick_flag = false;
        elseif Global.M == 10
            disable_flag = false;
            crowding_pick_flag = true;
        end
        normalization_str = 'normalize';
        reference_generation_mode = 'normal';
    end
    if strcmp(problem, 'MaF9')
        if norm(stable_threshold) < 1e-2
            stable_threshold = [5, 5, 5];
        end
        if Global.M == 5
            disable_flag = false;
            crowding_pick_flag = false;
        elseif Global.M == 10
            disable_flag = false;
            crowding_pick_flag = false;
        end
        normalization_str = '';
        reference_generation_mode = 'normal';
    end
    if strcmp(problem, 'MaF10')
        crowding_pick_flag = false;
        disable_flag = true;
        LARGE_ARCHIVE_FLAG = false;
    end
    if strcmp(problem, 'MaF11')
        crowding_pick_flag = false;
        disable_flag = true;
    end
    if strcmp(problem, 'MaF12')
        crowding_pick_flag = false;
        disable_flag = true;
    end
    if strcmp(problem, 'MaF13')
        crowding_pick_flag = false;
        disable_flag = false;
        reserved_evaluations = 100000;
        normalization_str = '';
    end
    if strcmp(problem, 'MaF14')%OK
        crowding_pick_flag = false;
        disable_flag = true;
        normalization_str = '';
    end
    if strcmp(problem, 'MaF15')
        if norm(stable_threshold) < 1e-2
            stable_threshold = [20, 20, 20];
        end
        disable_flag = false;
        crowding_pick_flag = false;
        normalization_str = '';
        if Global.M == 5
            reserved_evaluations = 300000;
        elseif Global.M == 10
            reserved_evaluations = 1850000;
        end
    end
    Z = generate_reference([], 0, Global, -inf, 'normal');
    reference_generation_intial_flag = true;
    Population = Global.Initialization(); Archive = Population; Sample = Population;
% end
end