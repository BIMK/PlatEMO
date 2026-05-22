function Population = ReMO_MapBack(Problem,LowDecs,MapCon,MapDiv,ConVars,DivVars)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if nargin < 4
        % Linear variable grouping
        ConVars = 1:round(Problem.D/2);
        DivVars = 1+round(Problem.D/2):Problem.D;
    end

    HighDecs = zeros(Problem.N,Problem.D);

    %% Mapping low-dimensional solution to high-dimensional space
    % Convegence
    ConDim  = size(MapCon,1);
    HighCon = LowDecs(:,1:ConDim)*MapCon;

    % Diversity
    DivDim  = size(MapDiv,1);
    HighDiv = LowDecs(:,1+ConDim:end)*MapDiv;

    
    HighDecs(:,ConVars) = HighCon;
    HighDecs(:,DivVars) = HighDiv;
    Lower = repmat(Problem.lower,Problem.N,1);
    Upper = repmat(Problem.upper,Problem.N,1);

    HighDecs   = min(max(HighDecs,Lower),Upper);
    Population = Problem.Evaluation(HighDecs);
end