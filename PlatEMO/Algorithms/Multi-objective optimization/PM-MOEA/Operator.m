function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,MaxP,MinP,Nonzero,PopDec)
% The operator of PM-MOEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,D]       = size(ParentMask);
    Parent1Mask = ParentMask(1:N/2,:);
    Parent2Mask = ParentMask(N/2+1:end,:);
    
    %% Delete illegal subspaces
    exist = false(size(MaxP,1),size(MinP,1));
    for i = 1 : size(MinP,1)
        exist(:,i) = all(MaxP(:,MinP(i,:)),2);
    end
    MaxP(~any(exist,2),:)  = [];
    MinP(~any(exist,1),:)  = [];

    %% Binary variation
    OffMask = false(N/2,D);
    if ~isempty(MaxP) && ~isempty(MinP)
        for i = 1 : N/2
            maxp = MaxP(randi(end),:);
            minp = MinP(randi(end),:);
            maxp = maxp & ~minp;
            OffMask(i,Nonzero(maxp)) = BinaryCrossover(Parent1Mask(i,Nonzero(maxp)),Parent2Mask(i,Nonzero(maxp)));
            OffMask(i,Nonzero(minp)) = true;
        end
    else
        OffMask = BinaryCrossover(Parent1Mask,Parent2Mask);
    end
    OffMask = BinaryMutation(OffMask);
    
    %% Real variation
    if any(Problem.encoding~=4)
        OffDec = OperatorGAhalf(Problem,ParentDec);
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(size(OffMask));
    end
    
    %% Remove duplicated solutions
    [~,uni] = unique(OffDec.*OffMask,'rows');
    OffDec  = OffDec(uni,:);
    OffMask = OffMask(uni,:);
    del            = ismember(OffDec.*OffMask,PopDec,'rows');
    OffDec(del,:)  = [];
    OffMask(del,:) = [];
end