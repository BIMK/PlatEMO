function lable = GetLable(Solution,non_dom)
% Get the label of each solution, 1 represents is

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    N = size(Solution,1);
    M = size(non_dom,1);
    Lables = zeros(N,M);

    %% Detect the dominance relation between each solution in Data and fns
    for i = 1 : N
        for j = 1 : M
            k = any(Solution(i,:)<non_dom(j,:)) - any(Solution(i,:)>non_dom(j,:));
            if k == 1 || k == 0
                Lables(i,j) = 1;
            end
        end
    end
    lable = zeros(1,N);
    lable(sum(Lables,2)==M) = 1;
end