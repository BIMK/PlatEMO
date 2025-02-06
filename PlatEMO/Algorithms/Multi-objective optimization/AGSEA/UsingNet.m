function [action, Fitness, Fitness3] = UsingNet(Problem,Fitness1,net,Memory,num_feature,max_act,Mask,action,Fitness3)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    output = zeros(1,max_act);
    if size(Memory,1) == 1      
    elseif rand < 0.3^(Problem.FE/Problem.maxFE) || size(Memory,1) < ceil(0.25*Problem.maxFE/100)
        action = randi([1 max_act],1);       
    else
        CurrentState = Memory(end,num_feature+2:end);
        for i = 1 : max_act
            input     = [CurrentState i]';
            output(i) = sim(net,input);
        end
        [~,action] = max(output);
    end     
    switch action
        case 1
            Fitness = Fitness1;
        case 2
            Fitness = sum(Mask==0);
            if size(Mask,1) == 1
                Fitness = Mask;
            end
        case 3
            Fitness = rand(1,Problem.D);
    end
end