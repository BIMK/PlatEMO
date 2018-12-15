function [KneePoints,Distance,r,t] = FindKneePoints(PopObj,FrontNo,MaxFNo,r,t,rate)
% Find all the knee points in each front

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(PopObj);

    %% Find the knee points in each front
    KneePoints = false(1,N);
    Distance   = zeros(1,N);
    for i = 1 : MaxFNo
        Current = find(FrontNo==i);
        if length(Current) <= M
            KneePoints(Current) = 1;
        else
            % Find the extreme points
            [~,Rank]   = sort(PopObj(Current,:),'descend');
            Extreme    = zeros(1,M);
            Extreme(1) = Rank(1,1);
            for j = 2 : length(Extreme)
                k = 1;
                Extreme(j) = Rank(k,j);
                while ismember(Extreme(j),Extreme(1:j-1))
                    k = k+1;
                    Extreme(j) = Rank(k,j);
                end
            end
            % Calculate the hyperplane
            Hyperplane = PopObj(Current(Extreme),:)\ones(length(Extreme),1);
            % Calculate the distance of each solution to the hyperplane
            Distance(Current) = -(PopObj(Current,:)*Hyperplane-1)./sqrt(sum(Hyperplane.^2));
            % Update the range of neighbourhood
            Fmax = max(PopObj(Current,:),[],1);
            Fmin = min(PopObj(Current,:),[],1);
            if t(i) == -1
                r(i) = 1;
            else
                r(i) = r(i)/exp((1-t(i)/rate)/M);
            end
            R = (Fmax-Fmin).*r(i);            
            % Select the knee points
            [~,Rank] = sort(Distance(Current),'descend');
            Choose   = zeros(1,length(Rank));
            Remain   = ones(1,length(Rank));
            for j = Rank
                if Remain(j)
                    for k = 1 : length(Current)     
                        if abs(PopObj(Current(j),:)-PopObj(Current(k),:)) <= R
                            Remain(k) = 0;
                        end
                    end
                    Choose(j) = 1;
                end
            end
            t(i) = sum(Choose)/length(Current);
            Choose(Rank(find(Choose(Rank)==1,1,'last'))) = 0;
            KneePoints(Current(Choose==1)) = 1;
        end
    end
end