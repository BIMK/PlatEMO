function subproblemList = ConstructLinearKRGOS(refPoints,Archive,nLinear,MaxObj,MinObj,W)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Get basic information
    % Get the number of dimension
    [DL,D] = size(refPoints);
    BU = ones(1,D);
    BD = zeros(1,D);
    [arcN,~] = size(Archive);
    % Check the number of reference point pairs
    if DL/2 ~= nLinear
        error('Check the number of reference points');
    end
    
    %% Construct linear space and build surrogate model
    % Tchebycheff approach with normalization
    Fitness = max(abs(Archive(:, D+1:end)-repmat(MinObj, arcN, 1))./repmat(MaxObj-MinObj,arcN,1).*W, [], 2);
    % Store linear space
    subproblemList = cell(1,nLinear);
    % Process
    for k = 0 : nLinear-1
        %% Construct linear space
        % Define subproblem struct
        subproblem = struct();
        % Get start point
        startPoint = refPoints(k*2+1,1:D);
        endPoint   = refPoints(k*2+2,1:D);
        subproblem.start = startPoint;
        % Calculate direction
        subproblem.direct = (endPoint-startPoint)./sqrt(sum((endPoint-startPoint).^2));
        % Find upper boundary and lower boundary of linear space
        [lub,llb]  = findBoundary(subproblem.start,subproblem.direct,BU,BD);
        subproblem.ub = lub;
        subproblem.lb = llb;
        
        %% Associate points to current linear space
        % Calculate the disance between start point and candidate points
        dist = pdist2(Archive(:,1:D),startPoint);
        % Get direction vector
        vector1 = repmat(subproblem.direct,arcN,1);
        % Get vector <point,startPoint>
        vector2 = Archive(:,1:D) - repmat(startPoint,arcN,1);
        % Calculate cos<vector1,vector2>
        MVL     = sum(vector2.^2,2);
        MVL(MVL==0) = inf;
        cosV    = sum(vector1.*vector2,2)./sqrt(MVL);
        % Calculate the distance between points and current line
        sinV    = sqrt(1-cosV.^2);
        allDist = sinV.*dist;
        % Calculate the transformed decision variables
        newX    = cosV.*dist;
        % associate process
        [sortDis,~] = sort(allDist,'ascend');
        trainSize   = 10;
        alpha  = sortDis(trainSize);
        select = allDist <= alpha;

        %% Build surrogate model
        % Radial basis function
        srgtOPTKRG = srgtsKRGSetOptions(newX(select),Fitness(select));
        srgtsKRG   = srgtsKRGFit(srgtOPTKRG);
        subproblem.fobj = @(Dec)srgtsKRGPredictor(Dec,srgtsKRG);

        subproblem.trainDec = newX(select);
        subproblem.maxTheta = mean(Fitness(select));
        subproblemList{k+1} = subproblem;
    end
end

function [linearBU,linearBD] = findBoundary(start,direct,BU,BD)
    lenMax = sqrt(sum((BU-BD).^2)) + 1e-8;
    
    %% Find the upper boundary of line
    ls = 0;
    le = lenMax;
    while true
        middle = (ls + le)/2;
        xNew   = start + middle*direct;
        if any(xNew>BU | xNew<BD) 
            le = middle;
        else
            ls = middle;
        end
        if le-ls < 1e-12
            break;
        end
    end
    linearBU = ls;
    
    %% Find the lower boundary of line
    ls = 0;
    le = -lenMax;
    while true
        middle = (ls + le)/2;
        xNew   = start + middle*direct;
        if any(xNew>BU | xNew<BD)
            le = middle;
        else
            ls = middle;
        end
        if ls-le < 1e-12
            break;
        end
    end
    linearBD = ls;
end