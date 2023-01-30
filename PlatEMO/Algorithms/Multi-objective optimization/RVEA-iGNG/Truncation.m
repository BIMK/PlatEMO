function Population = Truncation(Population,MaxSize)
% Limit the size of final popualtion in RVEA*

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu

    PopObj = Population.objs;
    zmax = max(PopObj,[],1);
    zmin = min(PopObj,[],1);

    PopObj= (PopObj - repmat(zmin,size(PopObj,1),1))./(zmax-zmin);
    M = size(PopObj,2);
    H = [eye(M-1)-ones(M-1)/M;-ones(1,M-1)/M];
    Pe = H*inv(H'*H)*H';
    f = PopObj*Pe';
    minf = min(f,[],1);
    maxf = max(f,[],1);
    f = (f-minf)./(maxf-minf);
    temp1 = f./sum(f,2);

    dis =  pdist2(temp1,temp1);
    dis(logical(eye(length(dis)))) = inf;
    Choose = true(1,length(Population));

    while sum(Choose) > MaxSize
        Remain   = find(Choose);
        Temp     = sort(dis(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Choose(Remain(Rank(1))) = false;
    end
    Population = Population(Choose);
end