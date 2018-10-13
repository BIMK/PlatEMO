function varargout = BT4(Operation,Global,input)
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
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            X     = input;
            [N,D] = size(X);
            
            I1         = 2:2:D;
            I2         = 3:2:D;
            Y          = X - sin(repmat(1:D,N,1)*pi/2/D);
            temp1      = X(:,1) < 0.25;
            temp2      = 0.25 <= X(:,1) & X(:,1) < 0.5;
            temp3      = 0.5 <= X(:,1) & X(:,1) < 0.75;
            temp4      = 0.75 <= X(:,1);
            X(temp1,1) = (1-(1-4*X(temp1,1)).^0.06)/4;
            X(temp2,1) = (1+(4*X(temp2,1)-1).^0.06)/4;
            X(temp3,1) = (3-(3-4*X(temp3,1)).^0.06)/4;
            X(temp4,1) = (3+(4*X(temp4,1)-3).^0.06)/4;
            PopObj(:,1) = X(:,1)         + sum(Y(:,I1).^2+(1-exp(-Y(:,I1).^2/1e-8))/5,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + sum(Y(:,I2).^2+(1-exp(-Y(:,I2).^2/1e-8))/5,2);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^0.5;
            varargout = {f};
    end
end