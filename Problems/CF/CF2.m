function varargout = CF2(Operation,Global,input)
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
            Global.lower    = [0,zeros(1,Global.D-1)-1];
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            X = input;
            D = size(X,2);
            
            J1 = 3:2:D;
            J2 = 2:2:D;
            PopObj(:,1) = X(:,1)         + 2*mean((X(:,J1)-sin(6*pi*repmat(X(:,1),1,length(J1))+repmat(J1,size(X,1),1)*pi/D)).^2,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean((X(:,J2)-cos(6*pi*repmat(X(:,1),1,length(J2))+repmat(J2,size(X,1),1)*pi/D)).^2,2);
            
            t = PopObj(:,2)+sqrt(PopObj(:,1))-sin(2*pi*(sqrt(PopObj(:,1))-PopObj(:,2)+1))-1;
            PopCon = -t./(1+exp(4*abs(t)));
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            X(:,1) = (0:1/(input-1):1)';
            X(:,2) = 1-sqrt(X(:,1));
            X(0<X(:,1) & X(:,1)<1/16 | 1/4<X(:,1) & X(:,1)<9/16,:) = [];
            varargout = {X};
    end
end