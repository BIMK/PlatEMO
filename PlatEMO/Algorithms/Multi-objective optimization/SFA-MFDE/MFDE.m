classdef MFDE
% Surrogate-assisted multifactorial differential evolution (MFDE) in SFA/MFDE

% This function was rewritten in the PlatEMO, building upon the MToP.
% We extend our sincere appreciation and respect to both parties involved.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
% Copyright (c) Yanchi Li. You are free to use the MToP for research
% purposes. All publications which use this platform should acknowledge
% the use of "MToP" or "MTO-Platform" and cite as "Y. Li, W. Gong, F. Ming,
% T. Zhang, S. Li, and Q. Gu, MToP: A MATLAB Optimization Platform for
% Evolutionary Multitasking, 2026, 6(3): Article No. 23, 1-30"
%--------------------------------------------------------------------------

% This function is written by Yuma Horaguchi

    properties (SetAccess = private)
        Problem;        % The PROBLEM class
        Task;           % Surrogate model
        RMP = 0.3;      % Random Mating Probability
        F   = 0.5;      % Scaling facter
        CR  = 0.9;      % Crossover rate
        T;              % The number of tasks
    end

    methods
        function obj = MFDE(Problem, Task, RMP, F, CR, T)
            obj.Problem = Problem;
            obj.Task    = Task;
            obj.RMP     = RMP;
            obj.F       = F;
            obj.CR      = CR;
            obj.T       = T;
        end

        function Population = run(obj, PopDec, omega)
            % Initialize population
            Population = obj.Initialization(PopDec);
            N          = length(Population) / obj.T;
            for gen = 1 : omega
                % Generate candidates
                [CandDec, CandMFFactor] = obj.Generation(Population);

                % Evaluation
                Candidate = obj.Evaluation(CandDec, CandMFFactor);

                % Update population
                Population = obj.SelectionMF([Population, Candidate], N);
            end
        end

        function Population = Initialization(obj, PopDec)
            % Calculate objectives
            PopObj = zeros(size(PopDec, 1), obj.T);
            for t = 1 : obj.T
                PopObj(:, t) = transpose(sim(obj.Task{t}, transpose(PopDec)));
            end

            % Calculate skill factor
            for TaskID = 1 : obj.T
                for SolID = 1 : size(PopDec, 1) / obj.T
                    Population(SolID + (TaskID - 1) * size(PopDec, 1) / obj.T) = SOLUTION(PopDec(SolID, :), PopObj(SolID, :), 0, TaskID);
                end
            end
        end

        function [OffDec, OffMFFactor] = Generation(obj, Population)
            OffDec      = [];
            OffMFFactor = Population.adds;

            for i = 1 : length(Population)
                Index = find(Population.adds == Population(i).add);
                [~, SrtID]  = sort(Population(Index).objs, 1);
                [~, BestID] = min(sum(SrtID(:, Population(i).add), 2));
                Index = setdiff(Index, i);
                x1    = Index(randi(length(Index)));
                if rand < obj.RMP
                    % Crossover with a solution having the different skill factor
                    Index = find(Population.adds ~= Population(i).add);
                    x2    = Index(randi(length(Index)));
                    OffMFFactor(i) = Population(x2).add;
                else
                    % Crossover with a solution having the same skill factor
                    Index = find(Population.adds == Population(i).add);
                    Index = setdiff(Index, [i, x1]);
                    x2    = Index(randi(length(Index)));
                    OffMFFactor(i) = Population(i).add;
                end
                % DE current-to-best/1
                NewDec  = Population(i).dec + obj.F * (Population(BestID).dec - Population(i).dec) + obj.F * (Population(x1).dec - Population(x2).dec);

                % Polynomial mutation
                Replace = rand(1, size(NewDec, 2)) > obj.CR;
                Replace(randi(size(NewDec, 2))) = false;
                ParDec  = Population(i).dec;
                NewDec(:, Replace) = ParDec(:, Replace);

                OffDec = [OffDec; min(max(NewDec, obj.Problem.lower), obj.Problem.upper)];
            end
        end

        function Offspring = Evaluation(obj, OffDec, OffMFFactor)
            % Calculate objectives
            OffObj = zeros(size(OffDec, 1), obj.T);
            for t = 1 : obj.T
                OffObj(:, t) = transpose(sim(obj.Task{t}, transpose(OffDec)));
            end

            % Convert to SOLUTION class
            for i = 1 : size(OffDec, 1)
                Offspring(i) = SOLUTION(OffDec(i, :), OffObj(i, :), 0, [OffMFFactor(i)]);
            end
        end

        function Next = SelectionMF(obj, Population, N)
            % Calculate facotrial ranks
            Next   = [];
            for t = 1 : obj.T
                Pop_t     = Population(Population.adds == t);
                Pop_tObj  = Pop_t.objs;
                [~, Rank] = sort(Pop_tObj(:, t));
                Next      = [Next, Pop_t(Rank(1 : N))];
            end
        end
    end
end