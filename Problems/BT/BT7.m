function varargout = BT7(Operation,Global,input)
% <problem> <BT>
% Biased Multiobjective Optimization and Decomposition Algorithm
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.M        = 2;
            Global.D        = 30;
            Global.lower    = [0,-ones(1,Global.D-1)];
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            X = input;
            D = size(X,2);
            
            I1 = 2:2:D;
            I2 = 3:2:D;
            Y  = X - repmat(sin(6*pi*X(:,1)),1,D);
            PopObj(:,1) = X(:,1)         + sum(Y(:,I1).^2+(1-exp(-Y(:,I1).^2/1e-3))/5,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + sum(Y(:,I2).^2+(1-exp(-Y(:,I2).^2/1e-3))/5,2);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^0.5;
            varargout = {f};
    end
end