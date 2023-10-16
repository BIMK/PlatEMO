function [Z, active_Z] = generate_reference(PopObj, original_active_number, Global, action, mode)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

global LEARNING_DISABLE_FLAG;
if strcmp(mode, 'strict')
    active_number = 0;
    PopObj = normalization(PopObj);
    active_Z = [];
    while active_number < original_active_number
        instructed_action = action;
        [Z, received_action] = generate_Z(Global, instructed_action, 'commit');
        Z = unique([Z; active_Z], 'rows');
        [~, Allocation] = pair(PopObj, Z, 'sin');
        active_index = unique(Allocation);
        active_Z = Z(active_index, :);
        active_number = size(active_Z, 1);
        fprintf('%d(+%d)\n', active_number, numel(active_index)); % TODO
        if received_action ~= instructed_action
            LEARNING_DISABLE_FLAG = true;
            warning('reference generation malfunctions');
            break;
        end
    end
elseif strcmp(mode, 'normal')
    Z = generate_Z(Global, action, 'commit');
    active_Z = [];
end
end

function varargout = generate_Z(Global, action, commit_str)
global MAX_REFERENCE_NUM_FLAG current_density reference_generation_intial_flag;
persistent incremental_sequence pointer;
if reference_generation_intial_flag == true
    MAX_REFERENCE_NUM_FLAG = false;
    incremental_sequence = density_sequence(Global.M);
    pointer = 1;
    reference_generation_intial_flag = false;
    % TRIM
    TRIM_INDEX = find(incremental_sequence < Global.N);
    incremental_sequence(TRIM_INDEX(1: end)) = [];
end
original_pointer = pointer;
pointer = pointer + action;
if pointer >= numel(incremental_sequence)
    MAX_REFERENCE_NUM_FLAG = true;
end
pointer = min(numel(incremental_sequence), pointer);
pointer = max(1, pointer);
real_action = pointer - original_pointer;

Z = UniformPoint(incremental_sequence(pointer), Global.M);
current_density = size(Z, 1);

if ~strcmp(commit_str, 'commit')
    pointer = original_pointer;
end

varargout{1} = Z;
if nargout >= 2
    varargout{2} = real_action;
end
end

function incremental_sequence = density_sequence(M)
incremental_sequence = [];
R_MAX = 1e6;
step = 1;
R = nchoosek(M + step - 1, step);
while R < R_MAX
    incremental_sequence = [incremental_sequence, R];
    step = step + 1;
    R = nchoosek(M + step - 1, step);
end
end