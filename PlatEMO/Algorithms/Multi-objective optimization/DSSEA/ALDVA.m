function rank_DV = ALDVA(Problem,Population,nSel,nPer)
% Detect the kind of each decision variable

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    [N,D] = size(Population.decs);
    ND    = NDSort(Population.objs,1) == 1;
    fmin  = min(Population(ND).objs,[],1);
    fmax  = max(Population(ND).objs,[],1);
    if any(fmax==fmin)
        fmax = ones(size(fmax));
        fmin = zeros(size(fmin));
    end
    
    %% Calculate the proper values of each decision variable
    Angle = zeros(D,nSel);
    RMSE  = zeros(D,nSel);
    
    for i = 1 : D
        drawnow('limitrate');
        Sample    = randi(N,1,nSel);    % Randomly select nSel solutions for perturbation
        % Generate several random solutions by perturbing the i-th dimension
        Decs      = repmat(Population(Sample).decs,nPer,1);
        Decs(:,i) = unifrnd(Problem.lower(i),Problem.upper(i),size(Decs,1),1);
        newPopu   = Problem.Evaluation(Decs);
        for j = 1 : nSel
            % Normalize the objective values of the current perturbed solutions
            Points      = newPopu(j:nSel:end).objs;
            After_sort  = sortrows(Points);
            Length(i,j) = norm(After_sort(1,:) - After_sort(end,:));%% Calculate the length after perturbation
            Points      = (Points-repmat(fmin,size(Points,1),1))./repmat(fmax-fmin,size(Points,1),1);
            Points      = Points - repmat(mean(Points,1),nPer,1);%%
            % Calculate the direction vector of the determining line
            [~,~,V] = svd(Points);
            Vector  = V(:,1)'./norm(V(:,1)');
            % Calculate the root mean square error
            error = zeros(1,nPer);
            for k = 1 : nPer
                error(k) = norm(Points(k,:)-sum(Points(k,:).*Vector)*Vector);
            end
            RMSE(i,j)  = sqrt(sum(error.^2));
            % Calculate the angle between the line and the hyperplane
            normal     = ones(1,size(Vector,2));
            sine       = abs(sum(Vector.*normal,2))./norm(Vector)./norm(normal);
            Angle(i,j) = real(asin(sine)/pi*180);%% asin: the larger it is, the more convergent it is
        end
    end

    Length = mean(Length,2);
    Angle  = mean(Angle,2);
    Length = (Length-min(Length))/(max(Length)-min(Length));
    Angle  = (Angle-min(Angle))/(max(Angle)-min(Angle));

    R           = Angle + Length;
    [~,rank_DV] = sort(R,'descend');
    [~,rank_DV] = sort(rank_DV);
    rank_DV     = rank_DV';
end