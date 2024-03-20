function Offspring = OperatorDE(Problem,Parent1,Parent2,Parent3,Parameter)
%OperatorDE - The operator of differential evolution.
%
%   Off = OperatorDE(Pro,P1,P2,P3) uses the operator of differential
%   evolution to generate offsprings for problem Pro based on parents P1,
%   P2, and P3. If P1, P2, and P3 are arrays of SOLUTION objects, then Off
%   is also an array of SOLUTION objects; while if P1, P2, and P3 are
%   matrices of decision variables, then Off is also a matrix of decision
%   variables, i.e., the offsprings are not evaluated. Each object or row
%   of P1, P2, and P3 is used to generate one offspring by P1 + 0.5*(P2-P3)
%   and polynomial mutation.
%
%	Off = OperatorDE(Pro,P1,P2,P3,{CR,F,proM,disM}) specifies the
%	parameters of operators, where CR and F are the parameters in
%	differental evolution, proM is the expectation of the number of mutated
%	variables, and disM is the distribution index of polynomial mutation.
%
%   Example:
%       Off = OperatorDE(Problem,Parent1,Parent2,Parent3)
%       Off = OperatorDE(Problem,Parent1.decs,Parent2.decs,Parent3.decs,{1,0.5,1,20})

%------------------------------- Reference --------------------------------
% H. Li and Q. Zhang, Multiobjective optimization problems with complicated
% Pareto sets, MOEA/D and NSGA-II, IEEE Transactions on Evolutionary
% Computation, 2009, 13(2): 284-302.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 4
        [CR,F,proM,disM] = deal(Parameter{:});
    else
        [CR,F,proM,disM] = deal(1,0.5,1,20);
    end
    if isa(Parent1(1),'SOLUTION')
        evaluated = true;
        Parent1   = Parent1.decs;
        Parent2   = Parent2.decs;
        Parent3   = Parent3.decs;
    else
        evaluated = false;
    end
    [N,D] = size(Parent1);

    %% Differental evolution
    Site = rand(N,D) < CR;
    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + F*(Parent2(Site)-Parent3(Site));

    %% Polynomial mutation
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    if evaluated
        Offspring = Problem.Evaluation(Offspring);
    end
end