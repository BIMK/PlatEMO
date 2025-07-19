function action = UsingNet(Problem,net,Memory,max_act)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    output = zeros(1,max_act);
    if rand < 0.3^(Problem.FE/Problem.maxFE) || size(Memory,1) < ceil(0.25*Problem.maxFE/100)
        action = randi([1 max_act],1);       
    else
        CurrentState = Memory(end,Problem.D + 1 +2:end);
        for i = 1 : max_act
            input     = [CurrentState i]';
            output(i) = sim(net,input);
        end
        [~,action] = max(output);
    end       
end