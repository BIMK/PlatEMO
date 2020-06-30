function x1 = FindSubspace(PopDec, epsilon)
% Pareto-optimal subspace learning in MaOEA-IT

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    x        = PopDec';
    avg      = mean(x,1);
    x        = x - repmat(avg,size(x,1),1);
    sigma    = x*x'/size(x,1);
    [~,S,~]  = svd(sigma);
    s        = diag(S);
    v_list   = cumsum(s)/sum(s);
    select_i = 1;
    while v_list(select_i) < epsilon
        select_i = select_i + 1;
    end
    x1   = PopDec;
    avg1 = round(mean(x1,1),2);
    x1(:,select_i+1:end) = repmat(avg1(select_i+1:end),size(x1,1),1);
end