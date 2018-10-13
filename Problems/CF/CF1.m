function varargout = CF1(Operation,Global,input)
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
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            X = input;
            D = size(X,2);
            
            J1 = 3:2:D;
            J2 = 2:2:D;
            PopObj(:,1) = X(:,1)   + 2*mean((X(:,J1)-repmat(X(:,1),1,length(J1)).^(0.5*(1+3*(repmat(J1,size(X,1),1)-2)/(D-2)))).^2,2);
            PopObj(:,2) = 1-X(:,1) + 2*mean((X(:,J2)-repmat(X(:,1),1,length(J2)).^(0.5*(1+3*(repmat(J2,size(X,1),1)-2)/(D-2)))).^2,2);
            
            PopCon = 1-PopObj(:,1)-PopObj(:,2)+abs(sin(10*pi*(PopObj(:,1)-PopObj(:,2)+1)));
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            PopObj(:,1) = (0:1/20:1)';
            PopObj(:,2) = 1-PopObj(:,1);
            varargout   = {PopObj};
    end
end