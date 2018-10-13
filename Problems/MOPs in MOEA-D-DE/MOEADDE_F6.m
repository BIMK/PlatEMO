function varargout = MOEADDE_F6(Operation,Global,input)
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
            Global.M        = 3;
            Global.M        = 3;
            Global.D        = 10;
            Global.lower    = [0,0,zeros(1,Global.D-2)-2];
            Global.upper    = [1,1,zeros(1,Global.D-2)+2];
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            X     = input;
            [N,D] = size(X);
            
            J1 = 4:3:D;
            J2 = 5:3:D;
            J3 = 3:3:D;
            PopObj(:,1) = cos(0.5*X(:,1)*pi).*cos(0.5*X(:,2)*pi) + 2*mean((X(:,J1)-2*repmat(X(:,2),1,length(J1)).*sin(repmat(2*pi*X(:,1),1,length(J1))+repmat(J1*pi/D,N,1))).^2,2);
            PopObj(:,2) = cos(0.5*X(:,1)*pi).*sin(0.5*X(:,2)*pi) + 2*mean((X(:,J2)-2*repmat(X(:,2),1,length(J2)).*sin(repmat(2*pi*X(:,1),1,length(J2))+repmat(J2*pi/D,N,1))).^2,2);
            PopObj(:,3) = sin(0.5*X(:,1)*pi)                     + 2*mean((X(:,J3)-2*repmat(X(:,2),1,length(J3)).*sin(repmat(2*pi*X(:,1),1,length(J3))+repmat(J3*pi/D,N,1))).^2,2);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,3);
            f = f./repmat(sqrt(sum(f.^2,2)),1,3);
            varargout = {f};
    end
end