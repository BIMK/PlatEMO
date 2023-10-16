function f = Distance_function(x, problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    upper = zeros(size(x)) + 10;
    lower = zeros(size(x));
    x = (upper - lower).* x + lower;
    if problem  ==1
        % Sphere function
        f = sum(x.^2,2);
    elseif problem == 2
        % Schwefel function
         f = max(abs(x),[],2);
    elseif problem == 3
        % Rastrigin function
        f = sum(x.^2-10.*cos(2.*pi.*x)+10,2);
    elseif problem == 4
        % Griewank function
        f = sum(x.^2,2)./4000 - prod(cos(x./repmat(sqrt(1:size(x,2)),size(x,1),1)),2) + 1;
    elseif problem == 5
        % Ackley function
    	f = 20-20.*exp(-0.2.*sqrt(sum(x.^2,2)./size(x,2)))-exp(sum(cos(2.*pi.*x),2)./size(x,2))+exp(1);
    end
    f(f<1e-8) = 0;
end