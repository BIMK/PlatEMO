function [Center,R] = adaptiveDivision(PopObj,K)
% Reference point adaption

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(PopObj);
    
    if K == 1
        Center = zeros(1,M);
        R      = inf;
        subNum = [1,N];
    else
       %% Detect the number of subregion
        % Calculate the distance between each solution
        fmin     = min(PopObj,[],1);
        fmax     = max(PopObj,[],1);
        PopObj   = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
        Distance = pdist2(PopObj,PopObj);    
        Distance(logical(eye(N))) = inf;
        radius   = max(min(Distance));
        % Detect subregion(s)
        Transformation = zeros(N,1);
        Remain         = find(Transformation==0);
        RegionID       = 1;
        while ~isempty(Remain)
            seeds      = find(~Transformation,1);
            Transformation(seeds) = RegionID;
            Remain     = find(Transformation==0);
            while true
                neighbors = sum(Distance(seeds,Remain)<=radius,1);
                seeds     = Remain(neighbors>=1);
                Transformation(seeds) = RegionID;
                Remain    = find(Transformation==0);
                if sum(neighbors)==0
                    break;
                end
            end
            RegionID = RegionID + 1;
        end
        %% Region division
        % Count the number of subregions of the true PF
        TrueNum = length(unique(Transformation));

        % Calculate the center point of each subregion
        Center = zeros(TrueNum,M);
        R      = ones(TrueNum,1);
        for i = 1 : TrueNum
            current     = Transformation==i;
            Center(i,:) = mean(PopObj(current,:));
            R(i)        = max(pdist2(PopObj(current,:),Center(i,:)));
        end

        % Select K points
        subNum  = tabulate(Transformation);
        subNum  = subNum(:,1:end-1);
        if TrueNum > K
            % Merging small subregions
            while sum(subNum(:,2)~=inf) > K
                [~,I] = min(subNum(:,2));
                Center(I,:) = inf(1,M);
                subNum(I,2) = inf;
                R(I)        = -inf;
                current     = find(Transformation == I);
                [~,T]       = min(pdist2(PopObj(current,:),Center),[],2);
                Transformation(current) = T;

                % Update reference point
                Idx = find(subNum(:,2)~=inf);
                for k =  1 : length(Idx)
                    Center(Idx(k),:) = mean(PopObj(Transformation == Idx(k),:));
                    R(Idx(k))        = max(pdist2(PopObj(Transformation == Idx(k),:),Center(Idx(k),:)))/sqrt(M-1);
                end
            end
        elseif TrueNum < K
            % Splite large subregions
            while sum(subNum(:,2)~=-inf) < K
                [~,I] = max(subNum(:,2));
                Center(I,:) = -inf(1,M);
                subNum(I,2) = -inf;
                R(I)        = -inf;
                current     = find(Transformation == I);
                [~,T1]      = max(pdist2(PopObj(current,:),PopObj(current(randi(length(current))),:)),[],1);
                [~,T2]      = max(pdist2(PopObj(current,:),PopObj(current(T1),:)),[],1);
                [~,T]       = min(pdist2(PopObj(current,:),PopObj(current([T1,T2],:),:)),[],2);
                ExistNum    = length(subNum(:,1));
                Transformation(current) = T + ExistNum;

                % Update reference point
                Center(ExistNum+1,:) = mean(PopObj(Transformation==ExistNum+1,:));
                Center(ExistNum+2,:) = mean(PopObj(Transformation==ExistNum+2,:));
                [R(ExistNum+1),R(ExistNum+2)] = deal(0.5*pdist2(Center(ExistNum+1,:),Center(ExistNum+2,:)));
                subNum(end+1,:) = [ExistNum+1,sum(T==1)];
                subNum(end+1,:) = [ExistNum+2,sum(T==2)];
            end
        end
    end
    
    % Select reference point
    select = abs(subNum(:,2)) ~= inf;
    Center = Center(select,:);
    R      = R(select);
end