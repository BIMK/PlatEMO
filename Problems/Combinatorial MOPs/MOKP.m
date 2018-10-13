function varargout = MOKP(Operation,Global,input)
% <problem> <Combinatorial MOP>
% Multi-objective knapsack problem
% operator --- EAbinary

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

persistent P W;

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 250;
            Global.operator = @EAbinary;
            
            [P,W] = Global.ParameterFile(sprintf('MOKP-M%d-D%d',Global.M,Global.D),randi([10,100],Global.M,Global.D),randi([10,100],Global.M,Global.D));
            
            PopDec    = randi([0,1],input,Global.D);
            PopDec    = MOKP('value',Global,PopDec);
            varargout = {PopDec};
        case 'value'
            X = input;
            C = sum(W,2)/2;
            [~,rank] = sort(max(P./W));
            for i = 1 : size(X,1)
                while any(W*X(i,:)'>C)
                    k = find(X(i,rank),1);
                    X(i,rank(k)) = 0;
                end
            end
            
            PopObj = repmat(sum(P,2)',size(X,1),1) - X*P';
            
            PopCon = [];
            
            varargout = {X,PopObj,PopCon};
        case 'PF'
            RefPoint  = sum(P,2)';
            varargout = {RefPoint};
        case 'draw'
            cla; Draw(input*P');
    end
end