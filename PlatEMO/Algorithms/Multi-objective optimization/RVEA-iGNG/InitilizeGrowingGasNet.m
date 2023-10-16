function net = InitilizeGrowingGasNet(V,Population,params)

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu

    N = params.N;
    MaxIt = params.MaxIt;
    L = params.L;
    epsilon_b = params.epsilon_b;
    epsilon_n = params.epsilon_n;
    alpha = params.alpha;
    delta = params.delta;
    T = params.T;

    PopObj = Population.objs;
    [NP,M]  = size(PopObj);
    PopObj = PopObj - repmat(min(PopObj,[],1),NP,1);
    Angle = acos(1-pdist2(PopObj,V,'cosine'));
    [~,associate] = min(Angle,[],2);
    valid = unique(associate);
    RefSize = size(valid,1);
    %% Initialization
    Ni = 2;
    w = zeros(Ni, M);
    if RefSize>=2
        for i = 1:Ni
            w(i,:) = V(valid(i,:),:);
        end
    else
        w(1:Ni,:) = V(randperm(N,Ni),:);
    end
    E = zeros(Ni,1);
    C = zeros(Ni, Ni);
    t = zeros(Ni, Ni);
    nx = 0;
    %% Loop
    for it = 1:MaxIt
        for kk = 3:RefSize
            % Select Input
            nx = nx + 1;
            x = V(valid(kk,:),:);

            % Competion and Ranking
            d = pdist2(x, w);
            [~, SortOrder] = sort(d);
            s1 = SortOrder(1);
            s2 = SortOrder(2);

            % Aging
            t(s1, :) = t(s1, :) + 1;
            t(:, s1) = t(:, s1) + 1;

            % Add Error
            E(s1) = E(s1) + d(s1)^2;

            % Adaptation
            w(s1,:) = w(s1,:) + epsilon_b*(x-w(s1,:));
            Ns1 = find(C(s1,:)==1);
            for j=Ns1
                w(j,:) = w(j,:) + epsilon_n*(x-w(j,:));
            end

            % Create Link
            C(s1,s2) = 1;
            C(s2,s1) = 1;
            t(s1,s2) = 0;
            t(s2,s1) = 0;

            % Remove Old Links
            C(t>T) = 0;
            nNeighbor = sum(C);
            AloneNodes = (nNeighbor==0);
            C(AloneNodes, :) = [];
            C(:, AloneNodes) = [];
            t(AloneNodes, :) = [];
            t(:, AloneNodes) = [];
            w(AloneNodes, :) = [];
            E(AloneNodes) = [];

            % Add New Nodes
            if mod(nx, L) == 0 && size(w,1) < N
                [~, q] = max(E);
                [~, f] = max(C(:,q).*E);
                r = size(w,1) + 1;
                w(r,:) = (w(q,:) + w(f,:))/2;
                C(q,f) = 0;
                C(f,q) = 0;
                C(q,r) = 1;
                C(r,q) = 1;
                C(r,f) = 1;
                C(f,r) = 1;
                t(r,:) = 0;
                t(:, r) = 0;
                E(q) = alpha*E(q);
                E(f) = alpha*E(f);
                E(r) = E(q);
            end

            % Decrease Errors
            E = delta*E;
        end
    %     PlotResults(w, C)
    end

    for ii = 1:size(w,1)
        ageSum(ii,:) = sum(t(ii,find(C(ii,:) == 1),:),2);
        ageSumBefore = ageSum;
        flag(ii,:) = 0;
    end
    net.w = w;
    net.E = E;
    net.C = C;
    net.t = t;
    net.nx = nx;
    net.ageSumBefore = ageSumBefore;
    net.flag = flag;
end