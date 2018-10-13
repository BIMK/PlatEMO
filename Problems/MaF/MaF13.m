function varargout = MaF13(Operation,Global,input)
% <problem> <MaF>
% A benchmark test suite for evolutionary many-objective optimization
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % This problem is P7
    switch Operation
        case 'init'
            Global.M          = 3;
            Global.D          = 5;
            Global.lower      = [zeros(1,2),zeros(1,Global.D-2)-2];
            Global.upper      = [ones(1,2),zeros(1,Global.D-2)+2];
            Global.operator   = @EAreal;
            Global.evaluation = max(1e5,1e4*Global.D);

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            X     = input;
            [N,D] = size(X);

            Y = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,N,1)*pi/D);
            PopObj(:,1) = sin(X(:,1)*pi/2)                   + 2*mean(Y(:,4:3:D).^2,2);
            PopObj(:,2) = cos(X(:,1)*pi/2).*sin(X(:,2)*pi/2) + 2*mean(Y(:,5:3:D).^2,2);
            PopObj(:,3) = cos(X(:,1)*pi/2).*cos(X(:,2)*pi/2) + 2*mean(Y(:,3:3:D).^2,2);
            PopObj(:,4:Global.M) = repmat(PopObj(:,1).^2+PopObj(:,2).^10+PopObj(:,3).^10+2*mean(Y(:,4:D).^2,2),1,Global.M-3);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,3);
            f = f./repmat(sqrt(sum(f.^2,2)),1,3);
            f = [f,repmat(f(:,1).^2+f(:,2).^10+f(:,3).^10,1,Global.M-3)];
            varargout = {f};
    end
end