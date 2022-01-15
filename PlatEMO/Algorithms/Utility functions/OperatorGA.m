function Offspring = OperatorGA(Parent,Parameter)
%OperatorGA - Crossover and mutation operators of genetic algorithm.
%
%   Off = OperatorGA(P) uses genetic operators to generate offsprings based
%   on the parents P. If P is an array of SOLUTION objects, then Off is
%   also an array of SOLUTION objects; while if P is a matrix of decision
%   variables, then Off is also a matrix of decision variables, i.e., the
%   offsprings are not evaluated. P is split into two subsets P1 and P2
%   with the same size, where each object or row of P1 and P2 is used to
%   generate two offsprings. Different operators are used for real, binary,
%   and permutation based encodings.
%
%   Off = OperatorGA(P,{proC,disC,proM,disM}) specifies the parameters of
%   operators, where proC is the probability of crossover, disC is the
%   distribution index of simulated binary crossover, proM is the
%   expectation of the number of mutated variables, and disM is the
%   distribution index of polynomial mutation.
%
%   Example:
%       Offspring = OperatorGA(Parent)
%       Offspring = OperatorGA(Parent.decs,{1,20,1,20})
%
%   See also OperatorGAhalf

%------------------------------- Reference --------------------------------
% [1] K. Deb, K. Sindhya, and T. Okabe, Self-adaptive simulated binary
% crossover for real-parameter optimization, Proceedings of the Annual
% Conference on Genetic and Evolutionary Computation, 2007, 1187-1194.
% [2] K. Deb and M. Goyal, A combined genetic adaptive search (GeneAS) for
% engineering design, Computer Science and informatics, 1996, 26: 30-45.
% [3] L. Davis, Applying adaptive algorithms to epistatic domains,
% Proceedings of the International Joint Conference on Artificial
% Intelligence, 1985, 162-164.
% [4] D. B. Fogel, An evolutionary approach to the traveling salesman
% problem, Biological Cybernetics, 1988, 60(2): 139-144.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 1
        [proC,disC,proM,disM] = deal(Parameter{:});
    else
        [proC,disC,proM,disM] = deal(1,20,1,20);
    end
    if isa(Parent(1),'SOLUTION')
        calObj = true;
        Parent = Parent.decs;
    else
        calObj = false;
    end
    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);
    Problem = PROBLEM.Current();
    
    switch Problem.encoding
        case 'binary'
            %% Genetic operators for binary encoding
            % Uniform crossover
            k = rand(N,D) < 0.5;
            k(repmat(rand(N,1)>proC,1,D)) = false;
            Offspring1    = Parent1;
            Offspring2    = Parent2;
            Offspring1(k) = Parent2(k);
            Offspring2(k) = Parent1(k);
            Offspring     = [Offspring1;Offspring2];
            % Bit-flip mutation
            Site = rand(2*N,D) < proM/D;
            Offspring(Site) = ~Offspring(Site);
        case 'permutation'
            %% Genetic operators for permutation based encoding
            % Order crossover
            Offspring = [Parent1;Parent2];
            k = randi(D,1,2*N);
            for i = 1 : N
                Offspring(i,k(i)+1:end)   = setdiff(Parent2(i,:),Parent1(i,1:k(i)),'stable');
                Offspring(i+N,k(i)+1:end) = setdiff(Parent1(i,:),Parent2(i,1:k(i)),'stable');
            end
            % Slight mutation
            k = randi(D,1,2*N);
            s = randi(D,1,2*N);
            for i = 1 : 2*N
                if s(i) < k(i)
                    Offspring(i,:) = Offspring(i,[1:s(i)-1,k(i),s(i):k(i)-1,k(i)+1:end]);
                elseif s(i) > k(i)
                    Offspring(i,:) = Offspring(i,[1:k(i)-1,k(i)+1:s(i)-1,k(i),s(i):end]);
                end
            end
        otherwise
            %% Genetic operators for real encoding
            % Simulated binary crossover
            beta = zeros(N,D);
            mu   = rand(N,D);
            beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
            beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
            beta = beta.*(-1).^randi([0,1],N,D);
            beta(rand(N,D)<0.5) = 1;
            beta(repmat(rand(N,1)>proC,1,D)) = 1;
            Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                         (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
            % Polynomial mutation
            Lower = repmat(Problem.lower,2*N,1);
            Upper = repmat(Problem.upper,2*N,1);
            Site  = rand(2*N,D) < proM/D;
            mu    = rand(2*N,D);
            temp  = Site & mu<=0.5;
            Offspring       = min(max(Offspring,Lower),Upper);
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                              (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
            temp = Site & mu>0.5; 
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                              (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    end
    if calObj
        Offspring = SOLUTION(Offspring);
    end
end