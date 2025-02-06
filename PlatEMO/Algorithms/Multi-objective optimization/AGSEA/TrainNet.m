function net = TrainNet(Problem,net,Memory,num_feature,max_act)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    inputn  = Memory(:,1:num_feature)';
    outputn = Memory(:,num_feature+1)';
    net.trainParam.epochs = 100;    
    net.trainParam.lr     = 0.001;         
    net.trainParam.goal   = 0.001;       
    net.trainParam.showWindow = 0;
    output = zeros(1,max_act);
    alpha  = 0.1;
    
    if size(Memory,1) == ceil(0.25*Problem.maxFE/100)
        net = newff(inputn,outputn,[10 10 10],{'tansig','purelin'},'trainlm'); 
    end
    if size(Memory,1) > ceil(0.25*Problem.maxFE/100) && mod(size(Memory,1),10) == 0
        MTrain = Memory(end-9:end,:);
        for i = 1 : size(MTrain,1)
            for j = 1 : max_act
                input = [MTrain(i,(num_feature+2):end) j]';
                output(j) = sim(net,input);
            end
            MTrain(i,num_feature+1) = MTrain(i,num_feature+1) + alpha*max(output);
        end
        inputn  = MTrain(:,1:num_feature)';
        outputn = MTrain(:,num_feature+1)';
        net     = train(net,inputn,outputn);
    end
end