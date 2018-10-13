function varargout = MOEADDE_F5(Operation,Global,input)
% <problem> <MOEADDE>
% Multiobjective Optimization Problems With Complicated Pareto Sets, MOEA/D
% and NSGA-II
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
            X     = input;
            [N,D] = size(X);
            
            J1   = 3:2:D;
            J2   = 2:2:D;
            X1J1 = repmat(X(:,1),1,length(J1));
            X1J2 = repmat(X(:,1),1,length(J2));
            PopObj(:,1) = X(:,1)         + 2*mean((X(:,J1)-(0.3*X1J1.^2.*cos(24*pi*X1J1+repmat(4*J1*pi/D,N,1))+0.6*X1J1).*cos(6*pi*X1J1+repmat(J1*pi/D,N,1))).^2,2);
            PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean((X(:,J2)-(0.3*X1J2.^2.*cos(24*pi*X1J2+repmat(4*J2*pi/D,N,1))+0.6*X1J2).*sin(6*pi*X1J2+repmat(J2*pi/D,N,1))).^2,2);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^0.5;
            varargout = {f};
    end
end