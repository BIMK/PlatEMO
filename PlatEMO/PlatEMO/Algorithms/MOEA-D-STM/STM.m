function Population = STM(Population,W,z,znad)
% Selection based on STM model

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N  = length(Population);
    NW = size(W,1);

    %% The modified Tchebycheff value of each solution on each subproblem
    g = zeros(N,NW);
    for i = 1 : N
        g(i,:) = max(repmat(abs(Population(i).obj-z),NW,1)./W,[],2)';
    end

    %% The perpendicular distance of each solution on each subproblem
    PopObj   = (Population.objs-repmat(z,N,1))./repmat(znad-z,N,1);
    Cosine   = 1 - pdist2(PopObj,W,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NW).*sqrt(1-Cosine.^2);
    
    %% STM selection
    Fp  = zeros(1,NW);
    FX  = zeros(1,N);
    Phi = false(NW,N);
    while any(Fp==0)
        RemainW  = find(Fp==0);
        i        = RemainW(randi(length(RemainW)));
        RemainX  = find(~Phi(i,:));
        [~,best] = min(g(RemainX,i));
        j        = RemainX(best);
        Phi(i,j) = true;
        if FX(j) == 0
            Fp(i) = j;
            FX(j) = i;
        elseif Distance(j,i) < Distance(j,FX(j))
            Fp(i)     = j;
            Fp(FX(j)) = 0;
            FX(j)     = i;
        end
    end
    Population = Population(Fp);
end