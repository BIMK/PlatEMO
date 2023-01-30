function Population = Mutation(Problem,Population)
% Mutation strategy of FROFI

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if ~any(all(Population.cons<=0,2))
        OffDec    = Population(randi(end)).dec;
        k         = randi(length(OffDec));
        OffDec(k) = unifrnd(Problem.lower(k),Problem.upper(k));
        Offspring = Problem.Evaluation(OffDec);
        [~,worst] = max(sum(max(0,Population.cons),2));
        if Population(worst).obj > Offspring.obj
            Population(worst) = Offspring;
        end
    end
end