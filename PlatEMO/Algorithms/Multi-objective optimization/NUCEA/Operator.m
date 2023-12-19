function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,Fitness,Mask)
% The operator of NUCEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [N,~]       = size(ParentDec);
    Parent1Dec  = ParentDec(1:floor(end/2),:);
    Parent2Dec  = ParentDec(floor(end/2)+1:floor(end/2)*2,:);
    Parent1Mask = ParentMask(1:floor(end/2),:);
    Parent2Mask = ParentMask(floor(end/2)+1:floor(end/2)*2,:);
    
    %% Clustering
    [~, index] = sort(Fitness);
    GroupSize  = ceil(mean(Mask,'all')*Problem.D);
    VaryGroup  = ones(1,Problem.D);
    start      = 1;
    GroupI     = 1;
    while true
        startend       = min(start+GroupI*GroupSize-1,Problem.D);
        ObjectiveIndex = start:startend;
        VaryGroup((index(ObjectiveIndex))) = GroupI;
        GroupI = GroupI + 1;
        start  = startend + 1;
        if start>Problem.D
            break;
        end
    end
    MaxGroup = max(VaryGroup);  

    %% Crossover and mutation for dec
    if any(Problem.encoding~=4)
        [OffDec,groupIndex,chosengroups] = GLP_OperatorGAhalf(Problem,Parent1Dec,Parent2Dec,4);	% 4 -- numberofgroups
        OffDec(:,Problem.encoding==4)    = 1;
    else
        OffDec = ones(size(Parent1Dec));
    end

    %% Crossover for mask
    OffMask = Parent1Mask;
    for i = 1 : N/2
        SelectedGroup = randi(MaxGroup,1);
        index = xor(Parent1Mask(i,:),Parent2Mask(i,:));     
        if rand < 0.5
            index = (SelectedGroup == VaryGroup) & index;           
            OffMask(i,index) = 0;
        else
            index = (SelectedGroup == VaryGroup) & index;           
            OffMask(i,index) = 1;
        end       
    end    
    if any(Problem.encoding~=4) && SelectedGroup < MaxGroup      
        chosenindex = groupIndex == chosengroups;
        for i = 1 : N/2            
            if rand < 0.5
                index = find(OffMask(i,:)&chosenindex(i,:));
                index = index(TS(-Fitness(index)));
                OffMask(i,index) = 0;
            else
                index = find(~OffMask(i,:)&chosenindex(i,:));
                index = index(TS(Fitness(index)));
                OffMask(i,index) = 1;
            end            
        end                    
    end       
end

function index = TS(Fitness)
% Binary tournament selection
    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(2,1,Fitness);
    end
end