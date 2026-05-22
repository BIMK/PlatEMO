function [Archive, ratio] = ExcuteAction(Problem, Archive, Action, Iter)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Form population
    switch Action
        case 1
            % Constrained optimization
            [Population1, ~, PopIndex] = FormPopulation(Archive, Action, Problem.N);
            Fitness1 = AdaFitness(Population1.objs, Population1.cons);
        case 2
            % Unconstrained optimization
            [~, Population2, PopIndex] = FormPopulation(Archive, Action, Problem.N);
            Fitness2 = AdaFitness(Population2.objs);
        case 3
            % Contrained and unconstrained optimization
            [Population1, Population2, PopIndex] = FormPopulation(Archive, Action, Problem.N);
            Fitness1 = AdaFitness(Population1.objs, Population1.cons);
            Fitness2 = AdaFitness(Population2.objs);
        otherwise
            error('Invalid action!');
    end
 
    %% optimization
    for i = 1 : Iter
        switch Action
            case 1
                % Constrained optimization
                Offspring1 = AdaSearch(Problem, Population1, Fitness1, false);
                [Population1,Fitness1] = AdaSelection([Population1,Offspring1],numel(Population1),true);
                if i == Iter
                    Population = Population1;
                end
            case 2
                % Unconstrained optimization
                Offspring2 = AdaSearch(Problem, Population2, Fitness2, false);
                [Population2,Fitness2] = AdaSelection([Population2,Offspring2],numel(Population2),false);
                if i == Iter
                    Population = Population2;
                end
            case 3
                % Contrained and unconstrained optimization
                Offspring1 = AdaSearch(Problem, Population1, Fitness1, true);
                Offspring2 = AdaSearch(Problem, Population2, Fitness2, true);
                [Population1,Fitness1] = AdaSelection([Population1,Offspring1, Offspring2],numel(Population1),true);
                [Population2,Fitness2] = AdaSelection([Population2,Offspring2, Offspring1],numel(Population2),false);
                if i == Iter
                    Population = [Population1, Population2];
                end            
            otherwise
                error('Invalid action!');
        end
    end
    if Action == 1 || Action == 2
        Archive(PopIndex) = Population;
    else
        [~, Index, ~] = unique(Population.decs, 'rows');
        Archive = Population(Index);
    end
    ratio = numel(Archive)/(1 + numel(Population));
end