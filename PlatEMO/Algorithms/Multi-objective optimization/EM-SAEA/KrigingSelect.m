function PopNew = KrigingSelect(PopDec,PopObj,V,mu,theta,PopCon)
% Kriging selection in K-RVEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);
    if nargin == 5
        CV = zeros(N,1);
    else
        CV = sum(max(0,PopCon),2);
    end

    [NVa,va] = NoActive(PopObj,V);
    NCluster = min(mu,size(V,1)-NVa);
    Va       = V(va,:);
    [IDX,~]  = kmeans(Va,NCluster);

    PopObj = PopObj - repmat(min(PopObj,[],1),size(PopObj,1),1);
    cosine = 1 - pdist2(Va,Va,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2);
    Angle  = acos(1-pdist2(PopObj,Va,'cosine'));
    [~,associate] = min(Angle,[],2);

    APD_S = ones(size(PopObj,1),1);
    for i = unique(associate)'
        current1 = find(associate==i);
        if ~isempty(current1)
            % Calculate the APD value of each solution
            APD = (1+size(PopObj,2)*theta*Angle(current1,i)/gamma(i)).*sqrt(sum(PopObj(current1,:).^2,2));
            % Select the one with the minimum APD value
            APD_S(current1,:) = APD;
        end
    end

    Cindex = IDX(associate); % Solution to cluster

    Next = zeros(NCluster,1);

    for i = unique(Cindex)'
        solution_Best1 = [];
        solution_Best2 = [];
        t = unique(associate(Cindex==i));
        % Calculate the APD value of each solution
        for j = 1 : size(t,1)
            currentS1 = find(associate==t(j)& CV==0);
            currentS2 = find(associate==t(j)& CV~=0);
            if ~isempty(currentS1)
                [~,id]         = min(APD_S(currentS1,:),[],1);
                solution_Best1 = [solution_Best1;currentS1(id)];
            elseif ~isempty(currentS2)
                [~,id]         = min(CV(currentS2));
                solution_Best2 = [solution_Best2;currentS2(id)];
            end
        end
        if ~isempty(solution_Best1)
            [~,best] = min(APD_S(solution_Best1,:),[],1);
            Next(i)  = solution_Best1(best);
        else
            [~,best] = min(CV(solution_Best2));
            Next(i)  = solution_Best2(best);
        end
    end
    index  = Next(Next~=0);
    PopNew = PopDec(index,:);
end