function PopNew = KrigingSelect(PopDec,PopObj,MSE,V,V0,NumV1,delta,mu,theta)
% Kriging selection in K-RVEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    [NumV2,~] = NoActive(PopObj,V0);
    [NVa,va]  = NoActive(PopObj,V);
    NCluster  = min(mu,size(V,1)-NVa);
    Va        = V(va,:);
    [IDX,~]   = kmeans(Va,NCluster);    

    PopObj = PopObj - repmat(min(PopObj,[],1),size(PopObj,1),1);
    cosine = 1 - pdist2(Va,Va,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2);
    Angle  = acos(1-pdist2(PopObj,Va,'cosine'));
    [~,associate] = min(Angle,[],2);
    Cindex = IDX(associate); % Solution to cluster

    Flag = NumV2-NumV1;
    Next = zeros(NCluster,1);

    for i = unique(Cindex)'
        current = find(Cindex==i);
        % Calculate the APD value of each solution
        if Flag<=delta
            APD = (1+size(PopObj,2)*theta*Angle(current,i)/gamma(i)).*sqrt(sum(PopObj(current,:).^2,2));
            % Select the one with the minimum APD value
            [~,best] = min(APD);
            Next(i)  = current(best);
        else
            Uncertainty = mean(MSE(current,:),2);
            [~,best]    = max(Uncertainty);
            Next(i)     = current(best);
        end
    end
    index  = Next(Next~=0);
    PopNew = PopDec(index,:);
end