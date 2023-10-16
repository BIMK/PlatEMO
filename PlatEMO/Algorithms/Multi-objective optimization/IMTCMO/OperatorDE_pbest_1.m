function Offspring = OperatorDE_pbest_1(Problem,Parent1,Parent2,Parent3,Parent4)
% The operator of differential evolution.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    %% Parameter setting
    [proM,disM] = deal(1,20);
    if isa(Parent1(1),'SOLUTION')
        evaluated  = true;
        Parent1 = Parent1.decs;
        Parent2 = Parent2.decs;
        Parent3 = Parent3.decs;
        Parent4 = Parent4.decs;
    else
        evaluated = false;
    end
    [N,D]   = size(Parent1);

    Fm = [0.6,0.8,1.0];
    CRm=[0.1,0.2,1.0];
    index = randi([1,length(Fm)],N,1);
    F = Fm(index);
    F = F';
    F = F(:, ones(1,D));
    index = randi([1,length(CRm)],N,1);
    CR = CRm(index);
    CR = CR';
    %% Differental evolution
    Site = rand(N,D) < repmat(CR,1,D);

    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + F(Site).*(Parent2(Site)-Offspring(Site) + Parent3(Site)-Parent4(Site));

    %% Polynomial mutation
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    if evaluated
        Offspring = Problem.Evaluation(Offspring);
    end
end