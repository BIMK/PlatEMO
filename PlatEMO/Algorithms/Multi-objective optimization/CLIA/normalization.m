function varargout = normalization(PopObj)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

global normalization_str
if strcmp(normalization_str, 'normalize')
    [N, ~] = size(PopObj);
    a = intercepts(PopObj);
    PopObj = (PopObj - repmat(min(PopObj, [], 1), size(PopObj, 1), 1)) ./ repmat(a,N,1);
end
varargout{1} = PopObj;
if nargout == 2
    varargout{2} = a;
end
end

function a = intercepts(PopObj)
[N, M] = size(PopObj);
%% Find the extreme points
[~, Choosed(1:M)] = min(PopObj, [], 1);
L2NormABO = zeros(N, M);
for i = 1 : M
    L2NormABO(:, i) = sum(PopObj(:, [1: i - 1, i + 1: M]) .^ 2, 2);
end
[~, Choosed(M + 1: 2 * M)] = min(L2NormABO, [], 1);
[~, Extreme] = max(PopObj(Choosed, :), [], 1);
Extreme = unique(Choosed(Extreme));
%% Calculate the intercepts
if length(Extreme) < M
    a = max(PopObj, [], 1);
else
    lastwarn('');
    Hyperplane = mldivide(PopObj(Extreme,:), ones(M, 1));
    [~, msgid] = lastwarn();
    if strcmp(msgid, 'MATLAB:nearlySingularMatrix')
        % error('error in normalize');
    end
    a = 1 ./ Hyperplane';
end
end