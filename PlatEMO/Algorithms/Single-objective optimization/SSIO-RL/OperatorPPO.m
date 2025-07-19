function Offspring = OperatorPPO(Population, index, Action, Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Parent  = Population(index).decs;
    Parent1 = Parent(1:floor(end/3),:);
    Parent2 = Parent(floor(end/3)+1:floor(end/3)*2,:); 
    Parent3 = Parent(floor(end/3)*2+1:end,:); 
    [N,D]   = size(Parent1);
    Lower   = repmat(Problem.lower,N,1);
    Upper   = repmat(Problem.upper,N,1);
    action  = Action(1:end,:);
    action  = min(max(action,-1),1);
    r1      = action(2:11);
    r2      = action(12:21);
    WeightP = action(22:end);
    WeightP = min(max(WeightP,1e-10),1);
    WeightProb = WeightP / sum(WeightP);
    Fit     = cumsum(WeightP);
    Fit     = Fit./max(Fit);
    type    = arrayfun(@(S)find(rand<=Fit,1),1:numel(Parent1));
    type    = reshape(type,size(Parent1));
    Qdec    = Parent1;
    for i = 1 : length(Fit)
        index = type == i;
        Qdec(index) = Parent3(index).*r1(i) + ...
                      Parent2(index).*r2(i) + ...
                      Parent1(index).*(1-r1(i)-r2(i));
    end
    Qdec       = min(max(Qdec,Lower),Upper);
    Site       = rand(N,D) < 1/D ;
    mu         = rand(N,D);
    temp       = Site & mu<=0.5;
    Qdec       = min(max(Qdec,Lower),Upper);
    Qdec(temp) = Qdec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*(1-(Qdec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^21).^(1/21)-1);
    temp       = Site & mu>0.5; 
    Qdec(temp) = Qdec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*(1-(Upper(temp)-Qdec(temp))./(Upper(temp)-Lower(temp))).^21).^(1/21));
    % Environmental selection
    Qdec = min(max(Qdec,Lower),Upper);
    Offspring = Problem.Evaluation(Qdec);
end

