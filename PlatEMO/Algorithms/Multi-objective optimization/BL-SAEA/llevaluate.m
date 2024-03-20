function [functionValue,equalityConstrVals,inequalityConstrVals]=llevaluate(llPop, Problem, ulMember)

    %This function evaluates the lower level objective values and constraints
    %for a set of lower level members corresponding to a given upper level member.
    noOfMembers = size(llPop,1);
    %Function call here

    equalityConstrVals = [];
    inequalityConstrVals = [];  
    llPopSize = (Problem.DL+1)*(Problem.DL+2)/2+Problem.DL;
    for i=1:noOfMembers
        if size(llPop,1) == size(ulMember,1) 
            llPopulation = Problem.EvaluationLower([ulMember(i,:), llPop(i,:)]);
            llPopCon = llPopulation.cons;
            llPopObj = llPopulation.objs;
            if sum(sum(isnan(llPopCon)))>0 || sum(sum(isnan(llPopObj)))>0
                llPopObj = llPopObj(~isnan(llPopObj));
                llPopCon = llPopCon(~isnan(llPopCon));
            end
            if(isempty(llPopCon))
                equalityConstrValsTemp = [];
                inequalityConstrValsTemp = [];
            else
                equalityConstrValsTemp = [];
                inequalityConstrValsTemp = llPopCon;
            end
            functionValue(i,:) = llPopObj;

        elseif size(ulMember,1) == 1 
            llPopulation = Problem.EvaluationLower([ulMember, llPop(i,:)]);
            llPopObj = llPopulation.objs;
            llPopCon = llPopulation.cons;
            if sum(sum(isnan(llPopCon)))>0 || sum(sum(isnan(llPopObj)))>0
                llPopObj = llPopObj(~isnan(llPopObj));
                llPopCon = llPopCon(~isnan(llPopCon));
            end
            if(isempty(llPopCon))
                equalityConstrValsTemp = [];
                inequalityConstrValsTemp = [];
            else
                equalityConstrValsTemp = [];
                inequalityConstrValsTemp = llPopCon;
            end
            functionValue(i,:) = llPopObj;
            
        else
            disp('Error in llTestProblem, size of llPop and ulMember mismatch');
        end
        
        if ~isempty(equalityConstrValsTemp)
            equalityConstrVals(i,:) = equalityConstrValsTemp;
        end
        if ~isempty(inequalityConstrValsTemp)
            inequalityConstrVals(i,:) = inequalityConstrValsTemp;
        end
    end
end