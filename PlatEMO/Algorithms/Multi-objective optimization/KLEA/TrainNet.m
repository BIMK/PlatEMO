function net = TrainNet(Problem,net,Memory,max_act)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    inputn  = Memory(:,1:(Problem.D + 1))';
    outputn = Memory(:,(Problem.D + 1)+1)';
    net.trainParam.epochs     = 100;    
    net.trainParam.lr         = 0.001;         
    net.trainParam.goal       = 0.001;       
    net.trainParam.showWindow = 0;
    output = zeros(1,max_act);
    alpha  = 1;
    
    if size(Memory,1) == ceil(0.25*Problem.maxFE/100)
        net = newff(inputn,outputn,[10 10 10],{'tansig','purelin'},'trainlm'); 
    end
    if size(Memory,1) > ceil(0.25*Problem.maxFE/100) && mod(size(Memory,1),10) == 0
        numRows = size(Memory, 1);      
        randIdx = randperm(numRows, 10); 
        MTrain  = Memory(randIdx, :); 
        for i = 1 : size(MTrain,1)
            for j = 1 : max_act
                input     = [MTrain(i,(Problem.D + 1+2):end) j]';
                output(j) = sim(net,input);
            end
            MTrain(i,Problem.D + 1 +1) = MTrain(i,Problem.D + 1 +1) + alpha*max(output);
        end
        inputn  = MTrain(:,1:(Problem.D + 1))';
        outputn = MTrain(:,Problem.D + 1 +1)';
        net     = train(net,inputn,outputn);
    end
end