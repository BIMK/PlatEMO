function [distance, Allocation] = pair(O, Z, type)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

BLOCK_SIZE = 10000;
blocks = floor(size(Z, 1) / BLOCK_SIZE);
B = mat2cell(Z, [BLOCK_SIZE * ones(1, blocks), size(Z, 1) - BLOCK_SIZE * blocks], size(Z, 2));
dA_set = {};
for i = 1: numel(B)
    z = B{i};
    dA = NaN(size(O, 1), 2);%[distance, Allocation]
    cor_index = (i - 1) * BLOCK_SIZE + 1: min(i * BLOCK_SIZE, size(Z, 1));
    Cosine = 1 - pdist2(O, z, 'cosine');
    if strcmp(type, 'distance')
        Error2Z = repmat(sqrt(sum(O .^ 2, 2)), 1, size(z, 1)) .* sqrt(1 - Cosine .^ 2);
    elseif strcmp(type, 'sin')
        Error2Z = sqrt(1 - Cosine .^ 2);
    end
    [dA(:, 1), dA(:, 2)] = min(Error2Z, [], 2);
    dA(:, 2) = cor_index(dA(:, 2));
    dA_set{i} = dA;
end
d = NaN(size(O, 1), numel(B)); A = NaN(size(O, 1), numel(B));
for i = 1: numel(B)
    dA = dA_set{i};
    d(:, i) = dA(:, 1); A(:, i) = dA(:, 2);
end
Allocation = NaN(size(O, 1), 1);
[distance, I] = min(d, [], 2);
for i = 1: size(d, 1)
    Allocation(i) = A(i, I(i));
end
end