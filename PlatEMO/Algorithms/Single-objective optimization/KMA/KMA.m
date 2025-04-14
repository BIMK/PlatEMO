classdef KMA < ALGORITHM
% <2022> <single> <real/integer> <large/none> <constrained/none>
% Komodo mlipir algorithm

%------------------------------- Reference --------------------------------
% S. Suyanto, A. A. Ariyanto, and A. F. Ariyanto. Komodo mlipir algorithm.
% Applied Soft Computing, 2022, 114: 108043.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            Big    = 1 : floor(0.5*Problem.N);
            Female = floor(0.5*Problem.N) + 1;
            Small  = floor(0.5*Problem.N)+2 : Problem.N;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                [~,rank]   = sort(FitnessSingle(Population));
                Population = Population(rank);
                % Movements of big males
                for i = 1 : length(Big)
                    OffDec1(i,:) = Population(Big(i)).dec + sum(repmat((-1).^((Population(Big).objs<Population(Big(i)).obj)|rand(length(Big),1)<0.5),1,Problem.D).*repmat(rand(length(Big),1),1,Problem.D).*(repmat(Population(Big(i)).dec,length(Big),1)-Population(Big).decs),1);
                end
                Pop             = [Population(Big),Problem.Evaluation(OffDec1)];
                [~,rank]        = sort(FitnessSingle(Pop));
                Population(Big) = Pop(rank(1:length(Big)));
                % Reproduction of female
                if rand < 0.5
                    OffDec2 = Population(Big(1)).dec + rand(1,Problem.D).*(Population(Female).dec-Population(Big(1)).dec);
                else
                    OffDec2 = Population(Female).dec + 0.1*(2*rand(1,Problem.D)-1).*(Problem.upper-Problem.lower);
                end
                Off = Problem.Evaluation(OffDec2);
                if Off.obj < Population(Female).obj
                    Population(Female) = Off;
                end
                % Movements of small males
                for i = 1 : length(Small)
                    OffDec3(i,:) = Population(Small(i)).dec + sum((rand(length(Big),Problem.D)<0.5).*rand(length(Big),Problem.D).*(Population(Big).decs-repmat(Population(Small(i)).dec,length(Big),1)),1);
                end
                Population(Small) = Problem.Evaluation(OffDec3);
            end
        end
    end
end