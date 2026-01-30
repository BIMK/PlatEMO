function models = TrainModel2(tr_x, tr_y, M, D)
% Train a RBFN Model by using the least square method to train weight.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yuanchao Liu (email:liuyuanchao@ise.neu.edu.cn)

    models = RBFMyself2(tr_x, tr_y, M, D);

    function models = RBFMyself2(tr_x, tr_y, M, D)
        ClusterNum = round(sqrt(M+D)+3);
        [models.Centers, models.Spreads, models.W2, models.B2] = clusterRBF(tr_x, tr_y,ClusterNum);
    end

    function [Centers, Spreads, W2, B2] = clusterRBF(SamIn, SamOut, ClusterNum)
        v       = size(SamIn,2);
        Overlap = 1.0;              % Overlap coefficient of hidden node
        SamNum  = size(SamIn,1);    % Number of all the samples
        nn      = 0;
        while nn==100 || nn==0
            index   = randi([1,SamNum],ClusterNum,1);
            Centers = SamIn(index,:);   %Initialise the clustering center
            n       = 1;
            while n < 100
                NumberInClusters = zeros(ClusterNum,1);         % Number of samples in each class, default is 0
                IndexInClusters  = zeros(ClusterNum,SamNum);    % Index of samples in each class
                % Classify all samples by the least distance principle
                for i = 1 : SamNum
                    AllDistance = dist(Centers,SamIn(i,:)');        % Calculate the distance between the i-th solution and each clustering center
                    [~,Pos]     = min(AllDistance);                 % Minimum distance,training input is the index of clustering center
                    NumberInClusters(Pos) = NumberInClusters(Pos) + 1;
                    IndexInClusters(Pos,NumberInClusters(Pos)) = i; % Stores the training indexes belonging to this class in turn
                end
                % Store the old clustering centers
                OldCenters = Centers;
                % Recalculate the clustering centers
                for i = 1 : ClusterNum
                    Index = IndexInClusters(i,1:NumberInClusters(i));   % Extract the training input index belonging to this class
                    Centers(i,:) = mean(SamIn(Index,:),1);              % Take the average of each class as the new clustering center
                end
                % Judge whether the old and new clustering centers are consistent
                EqualNum = sum(sum(Centers==OldCenters));   % Centers and Old Centers are subtracted from each other to sum up all corresponding bits
                if EqualNum == v*ClusterNum                 % The old and new clustering centers are consistent?
                    break;
                end
                n = n + 1;
            end
            nn = n;
        end
        % Calculate the spread constant (width) of each hidden node?
        AllDistances = dist(Centers,Centers');  % Calculate the distance between hidden node data centers (square matrix of ClusterNum dimension, symmetric matrix)?
        Maximum      = max(max(AllDistances));  % Find the largest distance?
        for i = 1 : ClusterNum                  % Replace the 0 on the diagonal with a larger value?
            AllDistances(i,i) = Maximum + 1;
        end
        Spreads = Overlap*min(AllDistances)';   % The minimum distance between hidden nodes is taken as the expansion constant.And convert it to a column vector?
        % Calculate the output weights of each hidden node
        Distance        = dist(Centers,SamIn');             % Calculate the distance between each sample input and each data center (Clusternum X Samnum matrix)
        SpreadsMat      = repmat(Spreads,1,SamNum);         % Clusternum X Samnum matrix
        HiddenUnitOut   = radbas(Distance./SpreadsMat);     % Calculate the hidden node output matrix;Radbas are radial basis transfer functions
        HiddenUnitOutEx = [HiddenUnitOut' ones(SamNum,1)]'; % Consider offsets (thresholds)
        W2Ex = SamOut'*pinv(HiddenUnitOutEx);               % Find the generalized output weight
        W2   = W2Ex(:,1:ClusterNum);                        % Output weight
        B2   = W2Ex(:,ClusterNum+1);                        % Offsets (thresholds)
    end
end