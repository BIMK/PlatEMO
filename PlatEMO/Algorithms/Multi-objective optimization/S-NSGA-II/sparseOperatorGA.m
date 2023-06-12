function Offspring = sparseOperatorGA(Problem, Parent, Parameter)
% Adapted from OperatorGA in PlatEMO by Ian Meyer Kropp

% This function is written by Ian Meyer Kropp

    %% Parameter setting
    if nargin > 1
        [proC,disC,proM,disM,proSM,disSM,mutation_method,crossover_method] = deal(Parameter{:});
    else
        [proC,disC,proM,disM,proSM,disSM,mutation_method,crossover_method] = deal(1,20,1,20,true,true);
    end
    
    calObj = false;

    if isa(Parent(1),'SOLUTION')
        calObj = true;
        Parent = Parent.decs;
    end

    % Check if any of the decision variables are non-real values 
    if any(ones(size(Problem.encoding)) ~= Problem.encoding)
        error('Only real encoding supported.');
    end

    Offspring = crossover_method(Parent, Problem.lower, Problem.upper, {proC,disC});

    mutation_params = {proM,disM, proSM,disSM};

    Offspring = mutation_method(Offspring, Problem.lower, Problem.upper, mutation_params);        

    if calObj
        Offspring = Problem.Evaluation(Offspring);
    end
end