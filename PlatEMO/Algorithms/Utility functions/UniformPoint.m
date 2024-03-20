function [W,N] = UniformPoint(N,M,method)
%UniformPoint - Generate a set of uniformly distributed points.
%
%   [W,L] = UniformPoint(N,M) returns approximately N uniformly distributed
%   points with M objectives on the unit hyperplane via the normal-boundary
%   intersection method with two layers. Note that the number of sampled
%   points L may be slightly smaller than the predefined size N due to the
%   need for uniformity.
%
%   [W,L] = UniformPoint(N,M,'ILD') returns approximately N uniformly
%   distributed points with M objectives on the unit hyperplane via the
%   incremental lattice design. Note that the number of sampled points L
%   may be slightly larger than the predefined size N due to the need for
%   uniformity.
%
%   W = UniformPoint(N,M,'MUD') returns exactly N uniformly distributed
%   points with M objectives on the unit hyperplane via the mixture uniform
%   design method.
%
%   [W,L] = UniformPoint(N,M,'grid') returns approximately N uniformly
%   distributed points with M objectives in the unit hypercube via the grid
%   sampling. Note that the number of sampled points L may be slighly
%   larger than the predefined size N due to the need for uniformity.
%
%   W = UniformPoint(N,M,'Latin') returns exactly N randomly distributed
%   points with M objectives in the unit hypercube via the Latin hypercube
%   sampling method.
%
%   Example:
%       [W,N] = UniformPoint(275,10)
%       [W,N] = UniformPoint(286,10,'ILD')
%       [W,N] = UniformPoint(102,10,'MUD')
%       [W,N] = UniformPoint(1000,3,'grid')
%       [W,N] = UniformPoint(103,10,'Latin')

%------------------------------- Reference --------------------------------
% [1] Y. Tian, X. Xiang, X. Zhang, R. Cheng, and Y. Jin, Sampling reference
% points on the Pareto fronts of benchmark multi-objective optimization
% problems, Proceedings of the IEEE Congress on Evolutionary Computation,
% 2018.
% [2] T. Takagi, K. Takadama, and H. Sato, Incremental lattice design
% of weight vector set, Proceedings of the Genetic and Evolutionary
% Computation Conference Companion, 2020, 1486-1494.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if nargin < 3
        method = 'NBI';
    end
    [W,N] = feval(method,N,M);
end

function [W,N] = NBI(N,M)
    H1 = 1;
    while nchoosek(H1+M,M-1) <= N
        H1 = H1 + 1;
    end
    W = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
    W = ([W,zeros(size(W,1),1)+H1]-[zeros(size(W,1),1),W])/H1;
    if H1 < M
        H2 = 0;
        while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
            H2 = H2 + 1;
        end
        if H2 > 0
            W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
            W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
            W  = [W;W2/2+1/(2*M)];
        end
    end
    W = max(W,1e-6);
    N = size(W,1);
end

function [W,N] = ILD(N,M)
    I = M * eye(M);
    W = zeros(1,M);
    edgeW = W;
    while size(W) < N
        edgeW = repmat(edgeW,M,1) + repelem(I,size(edgeW,1),1);
        edgeW = unique(edgeW,'rows');
        edgeW(min(edgeW,[],2)~=0,:) = [];
        W = [W+1;edgeW];
    end
    W = W./sum(W,2);
    W = max(W,1e-6);
    N = size(W,1);
end

function [W,N] = MUD(N,M)
    X = GoodLatticePoint(N,M-1).^(1./repmat(M-1:-1:1,N,1));
    X = max(X,1e-6);
    W = zeros(N,M);
    W(:,1:end-1) = (1-X).*cumprod(X,2)./X;
    W(:,end)     = prod(X,2);
end

function [W,N] = grid(N,M)
    gap = linspace(0,1,ceil(N^(1/M)));
    eval(sprintf('[%s]=ndgrid(gap);',sprintf('c%d,',1:M)))
    eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    N = size(W,1);
end

function [W,N] = Latin(N,M)
    [~,W] = sort(rand(N,M),1);
    W = (rand(N,M)+W-1)/N;
end

function Data = GoodLatticePoint(N,M)
    hm           = find(gcd(1:N,N)==1);
    udt          = mod((1:N)'*hm,N); 
    udt(udt==0)  = N;
    nCombination = nchoosek(length(hm),M);
    if nCombination < 1e4
        Combination = nchoosek(1:length(hm),M);
        CD2 = zeros(nCombination,1);
        for i = 1 : nCombination
            UT     = udt(:,Combination(i,:));
            CD2(i) = CalCD2(UT);
        end
        [~,minIndex] = min(CD2);
        Data = udt(:,Combination(minIndex,:));
    else
        CD2 = zeros(N,1);
        for i = 1 : N
            UT     = mod((1:N)'*i.^(0:M-1),N);
            CD2(i) = CalCD2(UT);
        end
        [~,minIndex] = min(CD2);
        Data = mod((1:N)'*minIndex.^(0:M-1),N);
        Data(Data==0) = N;
    end
    Data = (Data-1)/(N-1);
end

function CD2 = CalCD2(UT)
    [N,S] = size(UT);
    X     = (2*UT-1)/(2*N);
    CS1 = sum(prod(2+abs(X-1/2)-(X-1/2).^2,2));
    CS2 = zeros(N,1);
    for i = 1 : N    
        CS2(i) = sum(prod((1+1/2*abs(repmat(X(i,:),N,1)-1/2)+1/2*abs(X-1/2)-1/2*abs(repmat(X(i,:),N,1)-X)),2));
    end
    CS2 = sum(CS2);
    CD2 = (13/12)^S-2^(1-S)/N*CS1+1/(N^2)*CS2;
end