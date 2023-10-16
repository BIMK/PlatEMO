function Offspring = BinaryMutation(Offspring)
% Unbalanced binary mutation

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,D] = size(Offspring);
    MOne  = mean(Offspring,2);
    r     = min(min(1/D,2*MOne),2*(1-MOne));
    rate1 = repmat(r./2./MOne,1,D);
    rate0 = repmat(r./2./(1-MOne),1,D);
    rate  = zeros(N,D);
    rate(Offspring)  = rate1(Offspring);
    rate(~Offspring) = rate0(~Offspring);
    exchange            = rand(N,D) < rate;
    Offspring(exchange) = ~Offspring(exchange);
end