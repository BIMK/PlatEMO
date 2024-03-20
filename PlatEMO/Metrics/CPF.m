function score = CPF(Population,optimum)
% <max> <multi/many> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <multimodal/none> <sparse/none> <dynamic/none>
% Coverage over Pareto front

%------------------------------- Reference --------------------------------
% Y. Tian, R. Cheng, X. Zhang, M. Li, and Y. Jin, Diversity assessment of
% multi-objective evolutionary algorithms: Performance metric and benchmark
% problems, IEEE Computational Intelligence Magazine, 2019, 14(3): 61-74.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.best.objs;
    if size(PopObj,2) ~= size(optimum,2)
        score = nan;
    elseif size(optimum,1) > 1
        % Normalization
        fmin    = min(optimum,[],1);
        fmax    = max(optimum,[],1);
        PopObj  = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(fmax-fmin,size(PopObj,1),1);
        optimum = (optimum-repmat(fmin,size(optimum,1),1))./repmat(fmax-fmin,size(optimum,1),1);
        % Map to the Pareto front
        [~,Close] = min(pdist2(PopObj,optimum),[],2);
        PopObj    = optimum(Close,:);
        % Calculate the indicator value
        VPF   = Coverage(map(optimum,optimum),inf);
        V     = Coverage(map(PopObj,optimum),VPF/size(PopObj,1));
        score = V./VPF;
    else
        fmin   = min(PopObj,[],1);
        fmax   = max(PopObj,[],1);
        PopObj = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(fmax-fmin,size(PopObj,1),1);
        score  = Coverage(map(PopObj,PopObj),1/size(PopObj,1));
    end
end

function y = map(x,PF)
% Project the points in an (M-1)-d manifold to an (M-1)-d unit hypercube

    [N,M] = size(x);
    x  = x - repmat((sum(x,2)-1)/M,1,M);
    PF = PF - repmat((sum(PF,2)-1)/M,1,M);
    x  = x - repmat(min(PF,[],1),size(x,1),1);
    x  = x./repmat(sum(x,2),1,M);
    x  = max(1e-6,x);
    y  = zeros(N,M-1);
    for i = 1 : N
        c = ones(1,M);
        k = find(x(i,:)~=0,1);
        for j = k+1 : M
            temp     = x(i,j)/x(i,k)*prod(c(M-j+2:M-k));
            c(M-j+1) = 1/(temp+1);
        end
        y(i,:) = c(1:M-1);
    end
    y = y.^repmat(M-1:-1:1,N,1);
end

function V = Coverage(P,maxv)
% Calculate the hypervolume of each point's monopolized hypercube

    [N,M] = size(P);
    L = zeros(N,1);
    for x = 1 : N
        P1      = P;
        P1(x,:) = inf;
        L(x)    = min(max(abs(P1-repmat(P(x,:),N,1)),[],2));
    end
    L     = min(L,maxv.^(1/M));
    Lower = max(0,P-repmat(L/2,1,M));
    Upper = min(1,P+repmat(L/2,1,M));
    V     = sum(prod(Upper-Lower,2));
end