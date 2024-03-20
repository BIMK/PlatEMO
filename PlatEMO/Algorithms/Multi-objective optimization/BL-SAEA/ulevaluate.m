function [functionValue equalityConstrVals inequalityConstrVals]=ulevaluate(ulPop, llPop, Problem)

    %This function evaluates the upper level objective values and constraints
    %for a set of upper level members and their corresponding lower level members.
    %Function call here
    noOfMembers = size(ulPop,1);
    
    equalityConstrVals = [];
    inequalityConstrVals = [];
    l=0;
    for i=1:noOfMembers
        Population = Problem.Evaluation([ulPop(i, :), llPop(i, :)]);
        
        llPopulation = Problem.EvaluationLower([ulPop(i, :), llPop(i, :)]);
        llPopCon = llPopulation.cons;
        if sum(sum(isnan(llPopCon)))>0 || sum(sum(isnan(llPop)))>0
            llPopCon = llPopCon(~isnan(llPopCon));
        end
        l=length(llPopCon);

        PopCon = Population.cons;
        PopObj = Population.objs; 
        functionValue(i,:)=PopObj(:,1);
        inequalityConstrValsTemp = PopCon(:,1:l);
        equalityConstrValsTemp = [];

        if ~isempty(equalityConstrValsTemp)
            equalityConstrVals(i,:) = equalityConstrValsTemp;
        end
        if ~isempty(inequalityConstrValsTemp)
            inequalityConstrVals(i,:) = inequalityConstrValsTemp;
        end
    end
end