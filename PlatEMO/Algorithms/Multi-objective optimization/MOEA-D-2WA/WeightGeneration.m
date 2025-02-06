function [W1, N] = WeightGeneration(N, AdN, M, thetaM)
%Generate a set of uniform randomly distributed points on the objective and
%constraint violation space

%--------------------------------------------------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao

    SampleNum = 5000;
    W1        = UniformPoint(N, M);
    W1        = [W1 zeros(size(W1, 1), 1)];
    W2        = CalWeight1(thetaM, M+1, SampleNum);
    W1        = [W1; [zeros(1, M) 1]];
    while size(W1, 1) < N + AdN
        index = find_index_with_largest_angle (W1, W2);
        W1(size(W1, 1) + 1, :) = W2(index, :);
        W2(index, :) = [];
    end	
    
    W1 = max(W1, 1e-6);
    N  = size(W1, 1);
end

function theta = CalWeight1(thetaM, M, SampleNum)
    theta           = rand(SampleNum, M);
    theta           = theta./repmat(sum(theta, 2), 1, size(theta, 2));
    SumWeiFormer    = 1 - theta(:, M);
    SumWeiAfter     = 1 - theta(:, M).*thetaM;
    theta(:, 1:M-1) = SumWeiAfter./SumWeiFormer.*theta(:, 1:M-1);
    theta(:, M)     = theta(:, M).*thetaM;
end
    
function index = find_index_with_largest_angle (W1, W2)
    Cosine    = pdist2(W2, W1, 'cosine');  
    Temp      = sort(Cosine,2, 'descend');
    [~, Rank] = sortrows(Temp);
    index     = Rank(1);
end