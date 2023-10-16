function Population = EnvironmentalSelection(Population, N, Z, Zmin, initialE, epsn)
% The environmental selection of DCNSGA-III

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao
    
    [~, nCon] = size(Population.cons);
    ConVio    = max(0, Population.cons);
    if sum(sum(ConVio<=epsn, 2)==nCon) > N
        %% Selection among epsilon-feasible solutions
        tmp        = sum(ConVio<=epsn, 2)==nCon;
        Population = Population(1:end, tmp);
        CV         = sum(max(0,Population.cons)./initialE, 2)./nCon;
        %% Non-dominated sorting based on objectives and constraint violation
        [FrontNo, MaxFNo]    = NDSort([Population.objs, CV], N);
        Next                 = false(1, length(FrontNo));
        Next(FrontNo<MaxFNo) = true;
        %% Select solutions in the last front
        Last               = find(FrontNo==MaxFNo);
        pop1               = [Population(Next).objs];
        pop2               = [Population(Last).objs];
        Choose             = LastSelection(pop1, pop2, N-sum(Next), Z, Zmin);
        Next(Last(Choose)) = true;
        Population         = Population(Next);
    else
        %% Selection including epsilon-infeasible solutions
        CV         = sum(max(0, Population.cons)./initialE, 2)./nCon;
        [~, rank]  = sort(CV);
        Population = Population(rank(1:N));
    end
end

function Choose = LastSelection(PopObj1, PopObj2, K, Z, Zmin)
    %% Select part of the solutions in the last front
    PopObj = [PopObj1; PopObj2];
    [N, M] = size(PopObj);
    N1     = size(PopObj1, 1);
    N2     = size(PopObj2, 1);
    NZ     = size(Z, 1);
    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1, M);
    w       = zeros(M) + 1e-6 + eye(M);
    for i = 1 : M
        [~, Extreme(i)] = min(max(PopObj./repmat(w(i, :), N, 1), [], 2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M, 1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj, [], 1)';
    end
    % Normalization
    PopObj = (PopObj - repmat(Zmin, N, 1))./repmat(a' - Zmin, N, 1);
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj, Z, 'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2, 2)), 1, size(Z, 1)).*sqrt(1 - Cosine.^2);
    % Associate each solution with its nearest reference point
    [d, pi]  = min(Distance', [], 1);
    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = zeros(1, NZ);
    if N1 > 0
        temp = tabulate(pi(1:N1));
        rho(temp(:, 1)) = temp(:, 2);
    end
    %% Environmental selection
    Choose  = false(1, N2);
    Zchoose = true(1, NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1 + 1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~, s] = min(d(N1 + I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end