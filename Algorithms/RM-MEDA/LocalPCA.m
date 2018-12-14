function [Model,probability] = LocalPCA(PopDec,M,K)
% Partitioning the population by Local PCA algorithm

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is modified from the code in
% http://dces.essex.ac.uk/staff/zhang/IntrotoResearch/RegEDA.htm

    [N,D] = size(PopDec);
    Model = struct('mean',   num2cell(PopDec(1:K,:),2),...  % The mean of the model
                   'PI',     eye(D),...                     % The matrix PI
                   'eVector',[],...                         % The eigenvectors
                   'eValue', [],...                         % The eigenvalues
                   'a',      [],...                         % The lower bound of the projections
                   'b',      []);                           % The upper bound of the projections
    
    %% Modeling
    for iter = 1 : 50
        % Calculte the distance between each solution and its projection in
        % affine principal subspace of each cluster
        distance = zeros(N,K);
        for k = 1 : K
            distance(:,k) = sum((PopDec-repmat(Model(k).mean,N,1))*Model(k).PI.*(PopDec-repmat(Model(k).mean,N,1)),2);
        end
        % Partition
        [~,partition] = min(distance,[],2);
        % Update the model of each cluster
        updated = false(1,K);
        for k = 1 : K
            oldMean = Model(k).mean;
            current = partition == k;
            if sum(current) < 2
                if ~any(current)
                    current = randi(N);
                end
                Model(k).mean    = PopDec(current,:);
                Model(k).PI      = eye(D);
                Model(k).eVector = [];
                Model(k).eValue  = [];
            else
                Model(k).mean    = mean(PopDec(current,:),1);
                [eVector,eValue] = eig(cov(PopDec(current,:)-repmat(Model(k).mean,sum(current),1)));
                [eValue,rank]    = sort(diag(eValue),'descend');
                Model(k).eValue  = real(eValue);
                Model(k).eVector = real(eVector(:,rank));
                Model(k).PI      = Model(k).eVector(:,M:end)*Model(k).eVector(:,M:end)';
            end
            updated(k) = ~any(current) || sqrt(sum((oldMean-Model(k).mean).^2)) > 1e-5;
        end
        % Break if no change is made
        if ~any(updated)
            break;
        end
    end

	%% Calculate the smallest hyper-rectangle of each model
    for k = 1 : K
        if ~isempty(Model(k).eVector)
            hyperRectangle = (PopDec(partition==k,:)-repmat(Model(k).mean,sum(partition==k),1))*Model(k).eVector(:,1:M-1);
            Model(k).a     = min(hyperRectangle,[],1);
            Model(k).b     = max(hyperRectangle,[],1);
        else
            Model(k).a = zeros(1,M-1);
            Model(k).b = zeros(1,M-1);
        end
    end
    
    %% Calculate the probability of each cluster for reproduction
    % Calculate the volume of each cluster
    volume = prod(cat(1,Model.b)-cat(1,Model.a),2);
    % Calculate the cumulative probability of each cluster
    probability = cumsum(volume/sum(volume));
end