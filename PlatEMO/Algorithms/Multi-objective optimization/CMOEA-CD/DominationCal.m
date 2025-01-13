function NonDominated = DominationCal(Population, Add)
% The dominance relation calculation of CMOEA-CD
% add = 0 --- Pareto dominance
% add = 1 --- contraint-Pareto dominance

    PopObj = Population.objs;
    PopObj = roundn(PopObj, -10);
    [N, M] = size(PopObj);
    Dominated = false(1, N);
    PopCon = Population.cons;
    Cons = sum(max(0,PopCon),2);
    for i = 1: N-1
        err = PopObj(i,:) - PopObj;
        eq = zeros(N, 1);
        max_err = max(err, [],  2);
        min_err = min(err, [],  2);
        for j = i + 1: N
            for k = 1: M
                if err(j, k) ~= 0
                    break
                end
                if k == M
                    eq(j) = 1;
                end
            end
            if Add == 0
                if eq(j) == 1
                    Dominated(j) = true;
                elseif min_err(j) >= 0 
                    Dominated(i) = true;
                elseif max_err(j) <= 0
                    Dominated(j) = true;
                end
            else
                if eq(j) == 1
                    if Cons(i) <= Cons(j)
                        Dominated(j) = true;
                    else
                        Dominated(i) = true;
                    end
                elseif min_err(j) >= 0 
                    if Cons(j) <= 0 || Cons(j) <= Cons(i)
                        Dominated(i) = true;
                    end
                elseif max_err(j) <= 0
                    if Cons(i) <= 0 || Cons(i) <= Cons(j)
                        Dominated(j) = true;
                    end
                end
            end
        end
    end
    NonDominated = ~Dominated;
end
