function Offspring = FiveOps(Problem,op,x,x1,x2,x3,x4,x5)
% Four different DE operators

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [CR,F,proM,disM,K] = deal(1,0.5,1,20,0.5);
    D     = length(x.dec);
    Lower = Problem.lower;
    Upper = Problem.upper;

    %% Differental evolution
    switch op
        case 1
            % DE/rand/1
            v = x.dec + F*(x1.dec-x2.dec);
        case 2
            % DE/rand/2
            v = x.dec + F*(x1.dec-x2.dec) + F*(x3.dec-x4.dec);
        case 3
            % DE/current-to-rand/2
            v = x.dec + K*(x.dec-x1.dec) + F*(x2.dec-x3.dec) + F*(x4.dec-x5.dec);
        case 4
            % DE/current-to-rand/1
            v = x.dec + K*(x.dec-x1.dec) + F*(x2.dec-x3.dec);
        case 5
            uni = rand(1,D);
            v = x.dec + uni.*(Upper-Lower);

    end
    Offspring = x.dec;
    Site      = rand(1,D) < (CR+(op>2));
    Offspring(Site) = v(Site);

    %% Polynomial mutation
    Site  = rand(1,D) < proM/D;
    mu    = rand(1,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = Problem.Evaluation(Offspring);
end