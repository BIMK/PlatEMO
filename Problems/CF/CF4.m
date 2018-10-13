function varargout = CF4(Operation,Global,input)
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
            Global.M        = 2;
            Global.M        = 2;
            Global.D        = 10;
            Global.lower    = [0,zeros(1,Global.D-1)-2];
            Global.upper    = [1,zeros(1,Global.D-1)+2];
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            X = input;
            D = size(X,2);
            
            J1 = 3:2:D;
            J2 = 2:2:D;
            Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
            h          = Y.^2;
            temp       = Y(:,2)<3/2*(1-sqrt(1/2));
            h(temp,2)  = abs(Y(temp,2));
            h(~temp,2) = 0.125+(Y(~temp,2)-1).^2;
            PopObj(:,1) = X(:,1)   + sum(h(:,J1),2);
            PopObj(:,2) = 1-X(:,1) + sum(h(:,J2),2);
            
            t = X(:,2)-sin(6*pi*X(:,1)+2*pi/D)-0.5*X(:,1)+0.25;
            PopCon = -t./(1+exp(4*abs(t)));
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1) = (0:1/(input-1):1)';
            f(:,2) = 1-f(:,1);
            temp1  = 0.5<f(:,1) & f(:,1)<=0.75;
            temp2  = 0.75<f(:,1);
            f(temp1,2) = -0.5*f(temp1,1)+3/4;
            f(temp2,2) = 1-f(temp2,1)+0.125;
            varargout  = {f};
    end
end