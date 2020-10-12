function DAEA(Global)
% <algorithm> <D>
% Duplication Analysis Based Evolutionary Algorithm

%------------------------------- Reference --------------------------------
% H. Xu, B. Xue, and M. Zhang, A duplication analysis based evolutionary
% algorithm for bi-objective feature selection, IEEE Transactions on
% Evolutionary Computation, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Hang Xu

    %% Initialization
    Population = InitialPop(Global);

    %% Optimization
    while Global.NotTermination(Population)
        Offspring = NicVariation(Population);
        Population = EnvironmentalSelection([Population, Offspring], Global.N);
    end
end

function Pop = InitialPop(Global)
    N = Global.N;
    D = Global.D;
    T = min(D, N * 3);
    if T < D
        Pop = zeros(N, D);
        for i = 1 : N
            k = randperm(T, 1);
            j = randperm(D, k);
            Pop(i, j) = 1;
        end
        Pop = INDIVIDUAL(Pop);
    else
        Pop = Global.Initialization();
    end
end

function Offspring = NicVariation(Population)
    Objs = Population.objs;    
    Decs = Population.decs;
    [N, D] = size(Decs);

    % Selecting parents
    T = max(4, ceil(N * 0.2));
    normObjs = (mapminmax(Objs', 0, 1))';
    ED = pdist2(normObjs, normObjs, 'euclidean');
    ED(logical(eye(length(ED)))) = inf;
    [~, INic] = sort(ED, 2);
    INic = INic(:, 1 : T);
    IP_1 = (1 : N);
    IP_2 = zeros(1, N);
    for i = 1 : N
        if rand < 0.8 % local mating
            IP_2(i) = INic(i, randi(T, 1));
        else % global mating
            IG = (1 : N);
            IG(i) = [];
            IP_2(i) = IG(randi(N - 1, 1));
        end
    end
    Parent_1 = Decs(IP_1, :);
    Parent_2 = Decs(IP_2, :);
    Offspring = Parent_1;

    % do crossover
    for i = 1 : N
        k = find(xor(Parent_1(i, :), Parent_2(i, :)));
        t = length(k);
        if t > 1
            j = k(randperm(t, randi(t - 1, 1)));
            Offspring(i, j) = Parent_2(i, j);
        end
    end

    % do mutation
    for i = 1 : N
        if rand < 0.2
            j1 = find(Offspring(i, :));
            j0 = find(~Offspring(i, :));
            k1 = rand(1, length(j1)) < 1 / (length(j1) + 1);
            k0 = rand(1, length(j0)) < 1 / (length(j0) + 1);
            Offspring(i, j1(k1)) = false;
            Offspring(i, j0(k0)) = true;
        else
            k = rand(1, D) < 1 / D;
            Offspring(i, k) = ~Offspring(i, k);
        end
    end

    % get unique Offspring and individuals (function evaluated)
    Offspring = unique(Offspring, 'rows');
    Offspring = INDIVIDUAL(Offspring);
end

function Population = EnvironmentalSelection(Population, N)
    % Get unique individuals in decision space
    [~, U_Decs, ~] = unique(Population.decs, 'rows');
    UP = Population(U_Decs);
    Objs = UP.objs;
    Decs = UP.decs;

    if length(UP) > N 
        % Calculate solution difference in decision space
        SD = pdist2(Decs, Decs, 'cityblock');
        SD(logical(eye(length(SD)))) = inf;

        % remove some duplicated solutions in objective space
        [U_Objs, ~, I_Objs] = unique(Objs, 'rows');
        duplicated = [];
        D = size(Decs, 2);
        for i = 1 : size(U_Objs, 1)
            j = find(I_Objs == i);
            if length(j) > 1
                t = sum(Decs(j(1), :));
                d = min(SD(j, j), [], 2) / 2;
                p = d / t;
                r = find(p < 0.8 - 0.6 * (t - 1) / (D - 1));
                if ~isempty(r)
                    duplicated = [duplicated; j(r(randperm(length(r), length(r) - 1)))];
                end
            end
        end

        % reset population
        if length(UP) - length(duplicated) > N
            UP(duplicated) = [];
            Objs = UP.objs;
        end

        % nondominated sorting
        [Front, MaxF] = NDSort(Objs, N); 
        Selected = Front < MaxF;
        Candidate = Front == MaxF;

        % Calculate crowding distance
        CD = CrowdingDistance(Objs, Front, MaxF);

        % select last front
        while sum(Selected) < N
            S = Objs(Selected, 1);
            IC = find(Candidate);
            [~, ID] = sort(CD(IC), 'descend');
            IC = IC(ID);
            C = Objs(IC, 1);
            Div_Vert = zeros(1, length(C));
            for i = 1 : length(C)
                Div_Vert(i) = length(find(S == C(i)));
            end
            [~, IDiv_Vert] = sort(Div_Vert);
            IS = IC(IDiv_Vert(1));
            % reset Selected and Candidate
            Selected(IS) = true;
            Candidate(IS) = false;
        end
        Population = UP(Selected);
    else
        Population = [UP, Population(randperm(length(Population), (N - length(UP))))];
    end
end

function CrowdDis = CrowdingDistance(PopObj,FrontNo, MaxF)
    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Front = find(FrontNo==MaxF);
    Fmax  = max(PopObj(Front,:),[],1);
    Fmin  = min(PopObj(Front,:),[],1);
    for i = 1 : M
        [~,Rank] = sortrows(PopObj(Front,i));
        CrowdDis(Front(Rank(1)))   = inf;
        CrowdDis(Front(Rank(end))) = inf;
        for j = 2 : length(Front)-1
            CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
        end
    end
end