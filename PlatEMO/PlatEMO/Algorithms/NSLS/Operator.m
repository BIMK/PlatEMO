function Population = Operator(Population,Parameter)
% <operator> <real>
% The local search operator in NSLS
% mu    --- 0.5 --- The mean value of the Gaussian distribution
% delta --- 0.1 --- The standard deviation of the Gaussian distribution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 1
        [mu,delta] = deal(Parameter{:});
    else
        [mu,delta] = deal(0.5,0.1);
    end
    [N,D] = size(Population.decs);
    
    %% Local search
    for i = 1 : N
        for d = 1 : D
            c = mu + delta*randn;
            k = randi(N,1,2);
            w = repmat(Population(i).dec,2,1);
            w(1,d) = w(1,d) + c*(Population(k(1)).dec(d) - Population(k(2)).dec(d));
            w(2,d) = w(2,d) - c*(Population(k(1)).dec(d) - Population(k(2)).dec(d));
            w = INDIVIDUAL(w);
           	for j = 1 : 2
                k(j) = any(w(j).obj<Population(i).obj) - any(w(j).obj>Population(i).obj);
            end
            if k(1) == -1 && k(2) == -1
                continue;
            elseif k(1) > k(2)
                k = 1;
            elseif k(1) < k(2)
                k = 2;
            else
                k = randi(2);
            end
            Population(i) = w(k);
        end
    end
end