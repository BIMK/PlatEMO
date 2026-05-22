function Offspring = GESoperating(Problem, Parent, numberOfGroups, typeOfGroups, Model, A1)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

    if isa(Parent(1),'SOLUTION')
        evaluated = true;
        Parent    = Parent.decs;
    else
        evaluated = false;
    end
    [N,D]     = size(Parent);
    Parent1   = Parent(randperm(size(Parent, 1)), :);
    Parent2   = Parent(randperm(size(Parent, 1)), :);
    Offspring = Parent;
    C         = size(A1.cons,2);

    %% grouping
    [outIndexList,~] = CreateGroups(numberOfGroups,Parent,D,typeOfGroups);  % Create group
    chosengroups     = randi(numberOfGroups,size(outIndexList,1),1);        % An individual only selects a kind of grouping variables for optimization
    Site             = outIndexList == chosengroups; 

    %% GDE operators
    [CR,F,proM,disM] = deal(0.9,0.5,1,20);
    temp             = Site & rand(N,D) < CR;   % This includes crossover operations
    Offspring(temp)  = Parent(temp) + F*(Parent1(temp)-Parent2(temp));

    %% Linked Polynomial mutation
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    mu    = rand(N,1);
    mu    = repmat(mu,1,D);
    temp  = Site & mu<=0.5; 
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

    %% Evaluation
    if evaluated
        PopObj    = Model(Offspring)';
        PopDec    = Offspring;
        PopCon    = zeros(Problem.N,C);
        PopAdd    = zeros(Problem.N,1);
        Offspring = SOLUTION(PopDec,PopObj,PopCon,PopAdd);
    end
end