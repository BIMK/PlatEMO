function [OffDec,OffMask,action] = Operator(action,Problem,ParentDec,ParentMask,Fitness,Mask,Memory)
% The operator of KLEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
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
    
    %% Crossover and mutation for dec
    if any(Problem.encoding~=4)
        [OffDec,~,~] = GLP_OperatorGAhalf(Problem,Parent1Dec,Parent2Dec,4);
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(size(Parent1Dec));
    end

    %% Variable clustering
    [~, index] = sort(Fitness);
    GroupSize  = ceil(mean(Mask,'all')*Problem.D);  
    div        = Problem.D*length(unique(Memory(max(end-10,1):end,Problem.D + 1)))^2;
    ms         = 100*mean(Mask,'all');
    prob       = ms/div;
    switch action
        case 1
            VaryGroup = index;
            MaxGroup  = max(VaryGroup);
        case 2
            VaryGroup = ceil(index./GroupSize);
            MaxGroup  = max(VaryGroup);
        case 3
            VaryGroup = kmeans(Fitness',2)';
            MaxGroup  = max(VaryGroup);  
        otherwise
            error("Action invalid");
    end

    %% Crossover for mask
    OffMask = Parent1Mask;
    for i = 1 : N/2
        SelectedGroup = randi(MaxGroup,1);
        index = xor(Parent1Mask(i,:),Parent2Mask(i,:));
        if rand < 0.5
            index = (SelectedGroup == VaryGroup) & index & rand(1,Problem.D) < 1;           
            OffMask(i,index) = 0;
        else
            index = (SelectedGroup == VaryGroup) & index & rand(1,Problem.D) < 1;           
            OffMask(i,index) = 1;
        end       
    end

    for i = 1 : N/2    
        if rand < 0.5
            index = rand(1,Problem.D) < prob;           
            OffMask(i,index) = 0;
        else
            index = rand(1,Problem.D) < prob;           
            OffMask(i,index) = 1;
        end       
    end   
end