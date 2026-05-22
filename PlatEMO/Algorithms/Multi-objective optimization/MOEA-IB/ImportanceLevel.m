function [AlphaCon, AlphaDiv, ConVars, DivVars] = ImportanceLevel(Problem,Population,nSel,nPer)
% Detect the kind of each decision variable

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,D] = size(Population.decs);
    ND    = NDSort(Population.objs,1) == 1;
    fmin  = min(Population(ND).objs,[],1);
    fmax  = max(Population(ND).objs,[],1);
    if any(fmax==fmin)
        fmax = ones(size(fmax));
        fmin = zeros(size(fmin));
    end

    lower = min(Population.decs,[],1);
    upper = max(Population.decs,[],1);

    %% Avoid extremely small value ranges
    scalePop    = upper-lower;
    scalePro    = Problem.upper-Problem.lower;
    RepairLower = max(lower-0.1*scalePro,Problem.lower);
    RepairUpper = min(upper+0.1*scalePro,Problem.upper);
    Site        = scalePop<0.1*scalePro;
    lower(Site) = RepairLower(Site);
    upper(Site) = RepairUpper(Site);
    

    %% Calculate the proper values of each decision variable
    Diversity   = zeros(D,nSel);
    Convergence = zeros(D,nSel);
    RMSE   = zeros(D,nSel);
    Angle  = zeros(D,nSel);
    Sample = randi(N,1,nSel);

    for i = 1 : D
        drawnow('limitrate');
        % Generate several random solutions by perturbing the i-th dimension
        Decs = repmat(Population(Sample).decs,nPer,1);
        for j = 1 : nSel
            % Perturb
            Decs(j:nSel:end,i) = linspace(lower(i),upper(i),nPer);
            newPopu = Problem.Evaluation(Decs(j:nSel:end,:));

            % Normalize the objective values of the current perturbed solutions
            Points = newPopu.objs;
            Dist   = pdist2(Points,Points,'euclidean');

            Points = (Points-repmat(fmin,size(Points,1),1))./repmat(fmax-fmin,size(Points,1),1);
            Points = Points - repmat(mean(Points,1),nPer,1);
            % Calculate the direction vector of the determining line
            [~,~,V] = svd(Points);
            Vector  = V(:,1)'./norm(V(:,1)');

            % Calculate the root mean square error
            error = zeros(1,nPer);
            for k = 1 : nPer
                error(k) = norm(Points(k,:)-sum(Points(k,:).*Vector)*Vector);
            end
            RMSE(i,j) = sqrt(sum(error.^2));

            % Calculate the angle between the line and the hyperplane
            LDist  = max(Dist(:));

            normal = ones(1,size(Vector,2));
            cosine     = abs(sum(Vector.*normal,2))./norm(Vector)./norm(normal);
            Angle(i,j) = real(acos(cosine)/pi*180);


            Convergence(i,j) = cosine*LDist;

            sine = sqrt(1-cosine^2); 
            Diversity(i,j) = sine*LDist;
        end
    end


    %% Detect the kind of each decision variable
    VariableKind = (mean(RMSE,2)<1e-2)';
    result       = kmeans(Angle,2)';
    if any(result(VariableKind)==1) && any(result(VariableKind)==2)
        if mean(mean(Angle(result==1&VariableKind,:))) < mean(mean(Angle(result==2&VariableKind,:)))
            VariableKind = VariableKind & result==1;
        else
            VariableKind = VariableKind & result==2;
        end
    end

    %% Detect the kind of each decision variable
    Convergence = mean(Convergence,2);
    Diversity   = mean(Diversity,2);  
    ConVars     = find(VariableKind);
    DivVars     = find(~VariableKind);
    AlphaCon    = HarmonicMean(Convergence,VariableKind);
    AlphaDiv    = HarmonicMean(Diversity,~VariableKind);
end

function HMean = HarmonicMean(X,Props)
    if all(X(Props)==0)
        HMean = 0;
    else
        N     = sum(Props);
        X     = X/max(X);
        X     = max(X(Props),eps);
        HMean = N./sum(1./ X);
    end
end