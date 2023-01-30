function [V,net,genFlag] = TrainGrowingGasNet(V,temp1,net,scale,params,Problem,wholeObj,genFlag,zmin)

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu

    %% Parameters
    N = params.N;
    MaxIt = params.MaxIt;
    L = params.L;
    epsilon_b = params.epsilon_b;
    epsilon_n = params.epsilon_n;
    alpha = params.alpha;
    delta = params.delta;
    T = params.T;


    t = net.t;
    w = net.w;
    E = net.E;
    C = net.C;
    nx = net.nx;
    ageSumBefore = net.ageSumBefore;
    flag = net.flag;

    % select the corner solutions using GNG
    output=[];
    fm =[];
    C = net.C;
    w= net.w;
    for i = 1:size(w,1)
        neighbor = find(C(i,:)==1);

        ageSum(i,:) = sum(t(i,neighbor,:),2);
        if ageSum(i,:) == ageSumBefore(i,:)
            flag(i,:) = flag(i,:) + 1;
        end

    end

    cw = w(output,:);
    d = pdist2(cw,temp1);

    ageSumBefore = ageSum;
    maxN = 1.5;

    gen    = ceil(Problem.FE/Problem.N);
    maxgen = ceil(Problem.maxFE/Problem.N);    
    if gen <= round(0.9*maxgen)
        maxIter = 1;
        maxPZ = maxN;
        if size(w,1) == round(maxN*N)
            %      PlotResults(w, C)
            % plot(w(:,1),w(:,2),'k*')
            % hold on
            % plot(temp1(:,1),temp1(:,2),'rs')
            [~,rankFlag] = sort(flag,'descend');
            r = rankFlag(1:round(maxN*N)-N);
            %     [~,index] = max(flag,[],1);
            %     r = index;
            C(r, :) = [];
            C(:, r) = [];
            t(r, :) = [];
            t(:, r) = [];
            w(r, :) = [];
            E(r) = [];
            ageSumBefore(r,:) = [];
            flag(r,:) = [];
            flag = zeros(N,1);

        end
    else
        if size(w,1) < round(maxN*N) && isempty(genFlag)
            maxPZ = maxN;
            maxIter = 1;
        else
            maxPZ = 1;
            maxIter = 0;
        end
        if size(w,1) == round(maxN*N)
            maxPZ = 1;
            maxIter = 0;
            genFlag = gen;
        end
    end
    % plot(w(:,1),w(:,2),'k*')

    if isempty(genFlag)
        for iter = 1:maxIter
            for kk = 1:size(temp1,1)
                nx = nx + 1;
                x = temp1(kk,:);
                w = w./sum(w,2);
                d = pdist2(x, w);
                [~, SortOrder] = sort(d);
                s1 = SortOrder(1);
                s2 = SortOrder(2);

                % Aging: the age of all neighbours of s1 is increased by 1

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
                if ~isempty(find(AloneNodes == true))
                    %         AloneNodes
                end
                C(AloneNodes, :) = [];
                C(:, AloneNodes) = [];
                t(AloneNodes, :) = [];
                t(:, AloneNodes) = [];
                w(AloneNodes, :) = [];
                E(AloneNodes) = [];
                ageSumBefore(AloneNodes,:) = [];
                flag(AloneNodes,:) = [];

                % Add New Nodes

                if mod(nx, L) == 0 && size(w,1) < round(maxPZ*N)
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
                    ageSumBefore(r,:) = 0;
                    flag(r,:) = 0;
                end

                % Decrease Errors
                E = delta*E;
            end
        end


        w = w./sum(w,2);
        net.w = w;
        net.E = E;
        net.C = C;
        net.t = t;
        net.nx = nx;
        net.ageSumBefore = ageSumBefore;
        net.flag = flag;
    end

    if isempty(genFlag)
        output = [];
        for i = 1:size(w,1)
            neighbor = find(C(i,:)==1);
            [~,minInd] = min(w([i neighbor],:),[],1);
            [~,maxInd] = max(w([i neighbor],:),[],1);
            if ~isempty(find(minInd==1))
                output = [output i];
            elseif ~isempty(find(maxInd==1))
                output = [output i];
            end

        end

        for t = 1:size(output,2)
            s = output(:,t);
            x = mean(w(find(C(s,:)==1),:),1);
            w(s,:) = w(s,:) - 1*(x-w(s,:));
        end

        t = any(w<0,2);
        l = find(t == true);
        if ~isempty(l)
            for tt = 1:size(l,1)
                w(l(tt,1),any(w(l(tt,1),:)<0,1)) = 0;
            end
        end


        V = w;
        V = V.*repmat(scale,size(V,1),1);
    end





    if gen == genFlag
        N = size(wholeObj,1);
        wholeObj1 = (wholeObj - repmat(zmin,N,1));
        wholeObj = wholeObj1./scale;
        temp2 = wholeObj1./sum(wholeObj1,2);


        output = [];
        for i = 1:size(w,1)
            neighbor = find(C(i,:)==1);
            [~,minInd] = min(w([i neighbor],:),[],1);
            [~,maxInd] = max(w([i neighbor],:),[],1);
            if ~isempty(find(minInd==1))
                output = [output i];
            elseif ~isempty(find(maxInd==1))
                output = [output i];
            end

        end

        for t = 1:size(output,2)
            s = output(:,t);
            x = mean(w(find(C(s,:)==1),:),1);
            w(s,:) = w(s,:) - 1*(x-w(s,:));
        end

        t = any(w<0,2);
        l = find(t == true);
        if ~isempty(l)
            for tt = 1:size(l,1)
                w(l(tt,1),any(w(l(tt,1),:)<0,1)) = 0;
            end
        end

        V = [w.*repmat(scale,size(w,1),1);temp2];
    end
end