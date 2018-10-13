function varargout = CF8(Operation,Global,input)
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
            Global.lower    = [0,0,zeros(1,Global.D-2)-4];
            Global.upper    = [1,1,zeros(1,Global.D-2)+4];
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            X = input;
            D = size(X,2);
            
            J1 = 4:3:D;
            J2 = 5:3:D;
            J3 = 3:3:D;
            Y  = X-2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D) +repmat(1:D,size(X,1),1)*pi/D);
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean(Y(:,J1).^2,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean(Y(:,J2).^2,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean(Y(:,J3).^2,2);
            
            PopCon = 1-(PopObj(:,1).^2+PopObj(:,2).^2)./(1-PopObj(:,3).^2)+4*abs(sin(2*pi*((PopObj(:,1).^2-PopObj(:,2).^2)./(1-PopObj(:,3).^2)+1)));
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            input  = ceil(input/5)*5;
            f      = zeros(input,3);
            f(:,3) = repmat(sin((0:1/(input/5-1):1).*pi/2)',5,1);
            for i = 0 : 4
                f(i*input/5+1:(i+1)*input/5,1) = sqrt(i/4*(1-f(i*input/5+1:(i+1)*input/5,3).^2));
            end
            f(:,2) = sqrt(roundn(1-f(:,1).^2-f(:,3).^2,-10));
            varargout = {f};
    end
end