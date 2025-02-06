function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,REAL,Parameter)
% The Operator

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 4
        [proC,~,proM,~] = deal(Parameter{:});
    else
        [proC,~,proM,~] = deal(1,20,1,20);
    end
    Parent1Mask = ParentMask(1:floor(end/2),:);
    Parent2Mask = ParentMask(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]       = size(Parent1Mask);
    
    %% One point crossover and bitwise mutation for Mask
    k = repmat(1:D,N,1) > repmat(randi(D,N,1),1,D);
    k(repmat(rand(N,1)>proC,1,D)) = false;
    OffMask    = Parent1Mask;
    OffMask(k) = Parent2Mask(k);
    Site       = rand(N,D) < proM/D;
    OffMask(Site) = ~OffMask(Site);
    
    %% Crossover and mutation for dec
    if REAL
        OffDec = OperatorGAhalf(Problem,ParentDec);
    else
        OffDec = ones(N,D);
    end
end

