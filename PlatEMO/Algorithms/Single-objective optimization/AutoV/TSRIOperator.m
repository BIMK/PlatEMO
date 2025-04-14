function Offspring = TSRIOperator(Problem,Weight,Fit,Population)
% Offspring generation based on TSRI operator with given weights
% o = r1*(upper-lower) + r2*x2 + (1-r2)*x1, r1~N(0,w1^2), r2~N(w3,w2^2)

    %% Decision variables of parents
    Parent1 = Population(1:end/2).decs;
    Parent2 = Population(end/2+1:end).decs;
    [N,D]   = size(Parent1);
    UdL     = repmat(Problem.upper-Problem.lower,N,1);

    %% Offspring generation
    type = arrayfun(@(S)find(rand<=Fit,1),1:N*D);
    type = reshape(type,N,D);
    r1   = randn(N,D);
    r2   = randn(N,D);
    Odec = Parent1;
    for i = 1 : length(Fit)
        index = type == i;
        Odec(index) = UdL(index).*(r1(index)*Weight(i,1)) + ...
                      Parent2(index).*(r2(index)*Weight(i,2)+Weight(i,3)) + ...
                      Parent1(index).*(1-r2(index)*Weight(i,2)-Weight(i,3));
    end
    Offspring = Problem.Evaluation(Odec);
end