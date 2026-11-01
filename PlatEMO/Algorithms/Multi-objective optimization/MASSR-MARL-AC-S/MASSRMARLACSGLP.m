function Offspring = MASSRMARLACSGLP(Problem,Parent1,Parent2,numberOfGroups)
%MASSRMARLACSGLP - Apply grouped crossover and polynomial mutation.
%
%   Parent1 and Parent2 are decision matrices with one parent per row.
%   Construct simulated-binary offspring and mutate one value-ranked
%   variable group per offspring.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [proC,disC,disM] = deal(1,20,20);
    [N,D]            = size(Parent1);

    % Simulated binary crossover
    beta                             = zeros(N,D);
    mu                               = rand(N,D);
    beta(mu<=0.5)                    = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)                     = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta                             = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5)              = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring                        = (Parent1+Parent2)/2 + beta.*(Parent1-Parent2)/2;

    % Polynomial mutation in one selected variable group
    Lower           = repmat(Problem.lower,N,1);
    Upper           = repmat(Problem.upper,N,1);
    groupIndex      = CreateGroups(numberOfGroups,Offspring,D);
    chosenGroup     = randi(numberOfGroups,N,1);
    Site            = groupIndex == chosenGroup;
    mu              = rand(N,1);
    mu              = repmat(mu,1,D);
    temp            = Site & mu <= 0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp) + (Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).* ...
        (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp            = Site & mu > 0.5;
    Offspring(temp) = Offspring(temp) + (Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).* ...
        (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end

function outIndexArray = CreateGroups(numberOfGroups,xPrime,numberOfVariables)
%CreateGroups - Group variables by their offspring values.
%
%   Sort each offspring's variables and assign consecutive rank intervals to
%   numberOfGroups groups; the final group receives any remainder.

    outIndexArray = [];
    noOfSolutions = size(xPrime,1);
    for sol = 1 : noOfSolutions
        varsPerGroup = floor(numberOfVariables/numberOfGroups);
        vars         = xPrime(sol,:);
        [~,I]        = sort(vars);
        outIndexList = ones(1,numberOfVariables);
        for i = 1 : numberOfGroups-1
            outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
        end
        outIndexList(I(((numberOfGroups-1)*varsPerGroup)+1:end)) = numberOfGroups;
        outIndexArray                                            = [outIndexArray;outIndexList];
    end
end