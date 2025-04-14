classdef GCNMOEA < ALGORITHM
% <2025> <multi> <real/integer>
% Graph convolutional network based multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% P. Yan, Y. Tian, and Y. Liu. An indicator-based multi-objective
% evolutionary algorithm assisted by improved graph convolutional networks.
% Swarm and Evolutionary Computation, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            [W1,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population     = Problem.Initialization();
            Z              = min(Population.objs,[],1);
            epsilon_k      = 0;
            [~,FrontNo,d2] = EnvironmentalSelection(Population,W1,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                graph_data = Population.decs;
                maxConnectionsPerNode = 15;
                if rand < 0.3
                    R = corrcoef(graph_data'); 
                    adjacencyMatrix      = zeros(size(R));
                    [sortedR,sortedRIdx] = sort(R(:),'descend');
                else
                    R = corrcoef(graph_data'); 
                    adjacencyMatrix      = zeros(size(R));
                    [sortedR,sortedRIdx] = sort(R(:),'ascend');
                end
                for idx = 1 : numel(sortedR)
                    i = mod(sortedRIdx(idx)-1, size(R, 1))+1;
                    j = ceil(sortedRIdx(idx)/size(R, 1));
                    if sum(adjacencyMatrix(i,:)) < maxConnectionsPerNode && sum(adjacencyMatrix(:,j)) < maxConnectionsPerNode && i ~= j
                        adjacencyMatrix(i, j) = 1;
                        adjacencyMatrix(j, i) = 1;
                    end
                end
                G              = graph(adjacencyMatrix,'OmitSelfLoops');
                nodeProperties = array2table(graph_data);
                G.Nodes        = nodeProperties;
                if rand < 0.4 || Problem.FE< 0.5*Problem.maxFE
                    for num = 1 : 50
                        % red_population
                        red_node = randi([1 numnodes(G)]);
                        red_population{num} = red_node;
    
                        % green_population
                        green_node     = neighbors(G,red_node);
                        num_green_node = length(green_node);
                        if num_green_node > 8
                            green_node = datasample(green_node, 8, 'Replace', false);
                        end
                        green_population{num} = green_node;
                        W    = rand(size(graph_data,2),size(graph_data,2));
                        a    = rand(2*size(graph_data,2),1);
                        flag = 1;

                        % black_population
                        for i = 1 : length(green_population{num})
                            black_population = neighbors(G,green_population{num}(i));
                            black_population(black_population == red_population{num}) = [];
                            if length(black_population) > 10
                                black_population = datasample(black_population, 10, 'Replace', false);
                            end
                            if length(black_population) > 1
                                index      = [green_population{num}(i);black_population];
                                adjacency  = adjacencyMatrix(index,index);
                                adj_list   = mat2adj_list(adjacency);
                                max_clique = BronKerbosch(adj_list,[],1:length(black_population),[]);
                                max_clique = max_clique(2:end);
                                black_population1 = black_population(max_clique);
                                if length(black_population1) < 2
                                    if length(black_population) == 2
                                        black_population2 = datasample(black_population, 2, 'Replace', false);
                                    elseif length(black_population) == 3
                                        black_population2 = datasample(black_population, 3, 'Replace', false);
                                    else
                                        black_population2 = datasample(black_population, 4, 'Replace', false);
                                    end
                                    
                                    black_population2(black_population2 == black_population1) = [];
                                else
                                    black_population2 = [];
                                end
                                if isempty(black_population1)
                                    black_population1 = datasample(black_population, 4, 'Replace', false);
                                end
                                black_population = [black_population1;black_population2];
                                index  = [green_population{num}(i);black_population];
                                Parent = Population(index);
                                if length(Parent) == 2
                                    Offspring = OperatorGA(Problem,Parent,{1,10,1,10});
                                else
                                    adjacency = adjacencyMatrix(index,index);
                                    if flag == 1
                                        [W,a] = optimize_weights(Parent, adjacency, Problem);
                                        flag=2;
                                    end
                                    [~,Offspring] = self_attention(Parent,adjacency,W,a',Problem);
                                end
                                
                                Z = min(Z,min(Offspring.objs));
            
                                g_old  = max(abs(Parent.objs-Z).*W1(index,:),[],2);
                                g_new  = max(abs(Offspring.objs-Z).*W1(index,:),[],2);
                                cv_old = overall_cv(Parent.cons);
                                cv_new = overall_cv(Offspring.cons) .* ones(length(index),1);
                                ii     = index(find(((g_old >= g_new) & (((cv_old <= epsilon_k) & (cv_new <= epsilon_k)) | (cv_old == cv_new)) | (cv_new < cv_old) ), 2));
                                if ~isempty(ii)
                                    Population(ii) = Offspring(ismember(index,ii));
                                end
                            elseif isscalar(black_population)
                                Parent1    = Population(green_population{num}(i));
                                Parent2    = Population(black_population);
                                Parent     = [Parent1 Parent2];
                                Offspring  = OperatorGA(Problem,Parent);
                                Population = EnvironmentalSelection([Population,Offspring],W1,Problem.N);
                            end
                            [Population,FrontNo,d2] = EnvironmentalSelection(Population,W1,Problem.N);
                        end
                    end
                    
                    pg_ranks = centrality(G,'pagerank');
                    G.Nodes.PageRank = pg_ranks;
                    degree_ranks     = centrality(G,'degree');
                    bn_ranks         = centrality(G,'betweenness');
                    G.Nodes.Degree   = degree_ranks;
                    G.Nodes.Betweenness = bn_ranks;
                    G.Nodes.Weight = 0.5*pg_ranks + 0.5*bn_ranks;

                    for num = 1 : 50
                        Parent1  = Population(red_population{num});
                        Parents2 = table2array(G.Nodes(green_population{num},:));
                        len      = length(green_population{num});
                        if len == 1 || len == 2
                            [~,rowIndex] = max(Parents2(:,Problem.D+4));
                            Parent2      = Population(green_population{num}(rowIndex));
                            Parent       = [Parent1 Parent2];
                            Offspring    = OperatorGA(Problem,Parent,{1,10,1,10});
                            Population(red_population{num}) = Offspring(1);
                            Population(green_population{num}(rowIndex)) = Offspring(2);
                        elseif len > 0
                            [~,rowIndex] = sort(Parents2(:,Problem.D+4),'descend');
                            Parent2      = Population(green_population{num}(rowIndex(1)));
                            Parent3      = Population(green_population{num}(rowIndex(2)));
                            Offspring    = OperatorDE(Problem,Parent1,Parent2,Parent3);
                            Population(red_population{num}) = Offspring;
                        end
                        [Population,FrontNo,d2] = EnvironmentalSelection(Population,W1,Problem.N); 
                    end
                else
                    MatingPool = TournamentSelection(2,Problem.N,FrontNo,d2);
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                    Population = EnvironmentalSelection1([Population,Offspring],W1,Problem.N,Problem.M);
                    FrontNo    = NDSort(Population.objs,inf);
                    d2         = DensityEstimate(Population,W1);
                end                           
            end
        end
    end
end

function adj_list = mat2adj_list(mat)
    adj_list = cell(size(mat, 1), 1);
    for i = 1 : size(mat, 1)
        adj_list{i} = find(mat(i, :) == 1);
    end
end

function max_clique = BronKerbosch(adj_list, R, P, X)
    if isempty(P) && isempty(X)
        max_clique = R;
    else
        max_clique = [];
        for node = find_degree(P, adj_list)
            new_R  = union(R, node);
            new_P  = intersect(P, adj_list{node});
            new_X  = intersect(X, adj_list{node});
            clique = BronKerbosch(adj_list, new_R, new_P, new_X);
            if length(clique) > length(max_clique)
                max_clique = clique;
            end
            P = setdiff(P, node);
            X = union(X, node);
        end
    end
end

function max_clique = optimized_BronKerbosch(adj_list, R, P, X)
    if isempty(P) && isempty(X)
        max_clique = R;
    else
        max_clique = [];
        pivot      = find_pivot(P, X, adj_list);
        for node = setdiff(find_degree(P, adj_list), pivot)
            new_R = union(R, node);
            new_P = intersect(P, adj_list{node});
            new_X = intersect(X, adj_list{node});
            if length(new_P) > length(max_clique)
                clique = optimized_BronKerbosch(adj_list, new_R, new_P, new_X);
                if length(clique) > length(max_clique)
                    max_clique = clique;
                end
            end
            P = setdiff(P, node);
            X = union(X, node);
        end
    end
end

function pivot = find_pivot(P, X, adj_list)
    union_P_X     = union(P, X);
    max_neighbour = -1;
    pivot = P(1);
    for i = 1 : length(union_P_X)
        node_neighour = length(intersect(adj_list{union_P_X(i)},P));
        if max_neighbour < node_neighour
            max_neighbour = node_neighour;
            pivot = union_P_X(i);
        end
    end
end

function nodes = find_degree(P, adj_list)
    degrees = cellfun(@(node) length(intersect(P, adj_list{node})), num2cell(P));
    [~,ind] = sort(degrees, 'descend');
    nodes   = P(ind);
end

function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end

function [bestW, besta] = optimize_weights(input, adjacency, Problem)
    num_particles  = 30;
    num_iterations = 100;
    T_max = 100;
    T_min = 1e-3;
    cooling_rate = 0.99;

    T        = T_max;
    MAXFE    = Problem.FE;
    dimW     = size(input.decs, 2);
    dima     = dimW * 2;
    RefPoint = max(input.objs,[],1) + 0.1;

    W  = rand(num_particles, dimW, dimW);
    a  = rand(num_particles, dima);
    vW = rand(num_particles, dimW, dimW);
    va = rand(num_particles, dima);

    pbestW = W;
    pbesta = a;

    pbestHV = zeros(num_particles, 1);
    for i = 1 : num_particles
        [~,input1] = self_attention(input,adjacency,squeeze(W(i, :, :)),a(i, :)',Problem);
        pbestHV(i) = CalHV(input1.objs, RefPoint);
    end
    [global_best_HV,idx] = max(pbestHV);
    gbestW = W(idx, :, :);
    gbesta = a(idx, :);
    w_max  = 0.9;
    w_min  = 0.4;

    for iter = 1 : num_iterations
        w = w_max - (w_max - w_min) * (iter/num_iterations);
        for i = 1 : num_particles
            vW(i, :, :) = w * vW(i, :, :) + 3.0 * rand * (pbestW(i, :, :) - W(i, :, :)) + 3.0 * rand * (gbestW - W(i, :, :));
            va(i, :)    = w * va(i, :) + 3.0 * rand * (pbesta(i, :) - a(i, :)) + 3.0 * rand * (gbesta - a(i, :));
            W(i, :, :)  = W(i, :, :) + vW(i, :, :);
            a(i, :)     = a(i, :) + va(i, :);

            [~,input1] = self_attention(input,adjacency,squeeze(W(i, :, :)),a(i, :)',Problem);
            current_HV = CalHV(input1.objs, RefPoint);

            delta_HV = current_HV - pbestHV(i);
            if delta_HV > 0 || exp(delta_HV / T) > rand
                pbestW(i, :, :) = W(i, :, :);
                pbesta(i, :)    = a(i, :);
                pbestHV(i)      = current_HV;
            end

            if current_HV > global_best_HV
                global_best_HV = current_HV;
                gbestW = W(i, :, :);
                gbesta = a(i, :);
            end
        end
        T = cooling_rate * T;
        if T <= T_min
            break;
        end
    end
    bestW = squeeze(gbestW);
    besta = gbesta;
    Problem.FE = MAXFE;
end

function [H,input] = self_attention(input,adjacency,W,a,Problem)
    dec = input.decs;
    X   = dec;
    XW  = X * W;
    N   = size(XW, 1);
    e   = zeros(N, N);
    for i = 1 : N
        for j = 1 : N
            if adjacency(i, j) == 1
                e(i, j) = exp(tanh([XW(i, :), XW(j, :)] * a));
            else
                e(i, j) = 0;
            end
        end
    end

    alpha = zeros(N, N);
    for i = 1 : N
        alpha(i, :) = e(i, :) / sum(e(i, :));
    end

    H = zeros(N, size(W, 2));
    for i = 1 : N
        H(i, :) = sum(alpha(i, :)' .* XW, 1);
    end
    upper_bound = Problem.upper;
    lower_bound = Problem.lower;

    upper_bound_expanded = repmat(upper_bound, size(H, 1), 1);
    lower_bound_expanded = repmat(lower_bound, size(H, 1), 1);

    H = max(min(H, upper_bound_expanded), lower_bound_expanded);
    for i = 1 : length(input)
        input(i) = Problem.Evaluation(H(i,1:size(input.decs, 2)));
    end
end

