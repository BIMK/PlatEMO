function varargout = CF10(Operation,Global,input)
% <problem> <CF>
% Multiobjective optimization Test Instances for the CEC 2009 Special
% Session and Competition
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
            Global.D        = 10;
            Global.lower    = [0,0,zeros(1,Global.D-2)-2];
            Global.upper    = [1,1,zeros(1,Global.D-2)+2];
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            X = input;
            D = size(X,2);
            
            J1 = 4:3:D;
            J2 = 5:3:D;
            J3 = 3:3:D;
            Y  = X-2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(4*Y(:,J1).^2-cos(8*pi*Y(:,J1))+1,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(4*Y(:,J2).^2-cos(8*pi*Y(:,J2))+1,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(4*Y(:,J3).^2-cos(8*pi*Y(:,J3))+1,2);
            
            PopCon = 1-(PopObj(:,1).^2+PopObj(:,2).^2)./(1-PopObj(:,3).^2)+sin(2*pi*((PopObj(:,1).^2-PopObj(:,2).^2)./(1-PopObj(:,3).^2)+1));
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,3);
            f = f./repmat(sqrt(sum(f.^2,2)),1,3);
            f(1e-5<f(:,1) & f(:,1)<sqrt((1-f(:,3).^2)/4) | sqrt((1-f(:,3).^2)/2)<f(:,1) & f(:,1)<sqrt(3*(1-f(:,3).^2)/4),:) = [];
            varargout = {f};
    end
end