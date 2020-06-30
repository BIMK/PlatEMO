function [next] = Truncation(F, P, remaining)
% Truncation of MyO-DEMR, based on IGD

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

% get sizes
[N1, ~] = size(F);
[N2, ~] = size(P);

% normalize objectives
range = max(F)-min(F);
range(range==0) = 1;
F = (F - repmat(min(F), N1, 1))./repmat(range, N1, 1);

% build utopian front
R = P + min( min(sum(F,2)-1,0) );

% calculate distance matrix
dist = pdist2(R,F);

% select based on IGD
order = randperm(N2);
next = zeros(remaining,1);
for i = 1:remaining
    p = rem(i,N2); % point in ref set
    if p == 0
        p = N2;
    end
    [~, indexes] = sort( dist(order(p),:) );
    for idx = indexes
        if ~any(idx == next)
            next(i) = idx;
            break;
        end
    end
end

end