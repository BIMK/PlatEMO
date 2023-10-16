function Offspring = Operator(Problem,Population)
% Differential evolution in FROFI

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    PopDec = Population.decs;
    [N,D]  = size(PopDec);
    CR = [0.1 0.2 1]';
    CR = repmat(CR(randi(end,N,1)),1,D);
    F  = [0.6 0.8 1]';
    F  = repmat(F(randi(end,N,1)),1,D);
    
    %% Parents
    [~,P] = sort(rand(N),2);
    P1    = PopDec(P(:,1),:);
    P2    = PopDec(P(:,2),:);
    P3    = PopDec(P(:,3),:);
    [~,B] = min(Population.objs);
    PB    = repmat(PopDec(B,:),N,1);
    
    %% Offspring generation
    Rand   = rand(N,D);
    k1     = repmat(rand(N,1)<0.5,1,D);
    k2     = ~k1 & rand(N,D)<CR;
    OffDec = PopDec;
    OffDec(k1) = PopDec(k1) + Rand(k1).*(P1(k1)-PopDec(k1)) + F(k1).*(P2(k1)-P3(k1));
    OffDec(k2) = P1(k2) + Rand(k2).*(PB(k2)-P1(k2)) + F(k2).*(P2(k2)-P3(k2));
    Offspring  = Problem.Evaluation(OffDec);
end