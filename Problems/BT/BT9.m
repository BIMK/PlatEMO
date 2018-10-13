function varargout = BT9(Operation,Global,input)
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
            Global.M        = 3;
            Global.M        = 3;
            Global.D        = 30;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            X     = input;
            [N,D] = size(X);
            
            I1 = 3:3:D;
            I2 = 4:3:D;
            I3 = 5:3:D;
            Y  = X - sin(repmat(1:D,N,1)*pi/2/D);
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + sum(Y(:,I1).^2+(1-exp(-Y(:,I1).^2/1e-9))/5,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + sum(Y(:,I2).^2+(1-exp(-Y(:,I2).^2/1e-9))/5,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + sum(Y(:,I3).^2+(1-exp(-Y(:,I3).^2/1e-9))/5,2);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,3);
            f = f./repmat(sqrt(sum(f.^2,2)),1,3);
            varargout = {f};
    end
end