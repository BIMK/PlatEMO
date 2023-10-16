function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,Score)
% Operator of MP-MMEA

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------  
    
    %% Parameter setting
    N = size(ParentDec,1);
    Parent1Mask = ParentMask(1:N/2,:);
    Parent2Mask = ParentMask(N/2+1:end,:);
    
    %% Crossover for mask
    OffMask = Parent1Mask;
    for i = 1 : N/2    
        index1 = find(Parent1Mask(i,:)&~Parent2Mask(i,:));
        index2 = find(~Parent1Mask(i,:)&Parent2Mask(i,:));
        p1 = 1./(1+exp(-Score(index1)));
        p2 = 1./(1+exp(-Score(index2)));
        idx1 = index1(p1<rand(size(p1)));
        idx2 = index2(p2>rand(size(p2)));          
        OffMask(i,idx1) = 0;
        OffMask(i,idx2) = 1;
    end
    
    %% Mutation for mask
    for i = 1 : N/2    
        if rand < 0.5
            index = find(OffMask(i,:));
            index = index(TS(Score(index)));
            OffMask(i,index) = 0;
        else
            index = find(~OffMask(i,:));
            index = index(TS(-Score(index)));
            OffMask(i,index) = 1;
        end
    end
    
    %% Crossover and mutation for dec
    OffDec = OperatorGAhalf(Problem,ParentDec);
end

function index = TS(Fitness)
% Binary tournament selection

    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(2,1,Fitness);
    end
end