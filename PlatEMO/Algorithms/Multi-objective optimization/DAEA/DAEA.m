classdef DAEA < ALGORITHM
% <multi> <binary>
% Duplication analysis based evolutionary algorithm

%------------------------------- Reference --------------------------------
% H. Xu, B. Xue, and M. Zhang, A duplication analysis based evolutionary
% algorithm for bi-objective feature selection, IEEE Transactions on
% Evolutionary Computation, 2021, 25(2): 205-218.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Hang Xu

    methods
        function main(Algorithm,Problem)
            %% Initialization
            Population = InitialPop(Problem);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                Offspring  = NicVariation(Problem,Population);
                Population = EnvironmentalSelection([Population, Offspring],Problem.N);
            end
        end
    end
end

function Pop = InitialPop(Problem)
    N = Problem.N;
    D = Problem.D;
    T = min(D, N * 3);
    if T < D
        Pop = false(N, D);
        for i = 1 : N
            k = randperm(T, 1);
            j = randperm(D, k);
            Pop(i, j) = 1;
        end
        Pop = Problem.Evaluation(Pop);
    else
        Pop = Problem.Initialization();
    end
end

function Offspring = NicVariation(Problem,Population)
    Objs  = Population.objs;    
    Decs  = Population.decs;
    [N,D] = size(Decs);

    % Selecting parents
    T = max(4, ceil(N * 0.2));
    normObjs = (Objs-repmat(min(Objs,[],1),N,1))./repmat(max(Objs,[],1)-min(Objs,[],1),N,1);
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
    Parent_1  = Decs(IP_1, :);
    Parent_2  = Decs(IP_2, :);
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

    % get unique offspring and individuals (function evaluated)
    Offspring = unique(Offspring, 'rows');
    Offspring = Problem.Evaluation(Offspring);
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
        CD = CrowdingDistance(Objs,Front);

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