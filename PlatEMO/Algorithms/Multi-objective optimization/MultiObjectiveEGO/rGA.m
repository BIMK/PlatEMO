function PopDec=rGA(fhandle,Problem)
% real-parameter genetic algorithm to find minimal objective value

%--------------------------------------------------------------------------
% This function is written by Youwei He (email: 1554748356@qq.com)

% the number of design variables
D = Problem.D;
% population size of GA
GAPopulationSize = 10*D;
% Generation of GA
GAGeneration = 100;
obj_max=Inf;
% the first GA generation, randomly generated
Offspring = repmat(Problem.upper-Problem.lower,GAPopulationSize,1).*...
            UniformPoint(GAPopulationSize,D,'Latin')+repmat(Problem.lower,GAPopulationSize,1);
% the GA process for optimizing the objective function
for gen = 1 :  GAGeneration
    obj_Offspring = feval(fhandle, Offspring);
    [~,index] = sort(obj_Offspring,'ascend');
    if obj_Offspring(index(1)) < obj_max
        Best = Offspring(index(1),:);
        obj_max   = obj_Offspring(index(1));
    end
    Parent    = Offspring(index(1:ceil(GAPopulationSize/2)),:);
    Offspring = [OperatorGA(Problem,Parent(TournamentSelection(2,size(Parent,1), ...
        obj_Offspring(index(1:ceil(GAPopulationSize/2)))),:));OperatorGA(Problem,Parent,{0.9,2,1/D,20})];
end
PopDec = Best;
end