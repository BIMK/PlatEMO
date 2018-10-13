function varargout = MOEADM2M_F6(Operation,Global,input)
% <problem> <MOEADM2M>
% Decomposition of a Multiobjective Optimization Problem into a Number of
% Simple Multiobjective Subproblems
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
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            X = input;
            
            t = X(:,3:end) - repmat(X(:,1).*X(:,2),1,size(X,2)-2);
            g = 2*sin(pi*X(:,1)).*sum(-0.9*t.^2+abs(t).^0.6,2);
            PopObj(:,1) = (1+g).*X(:,1).*X(:,2);
            PopObj(:,2) = (1+g).*X(:,1).*(1-X(:,2));
            PopObj(:,3) = (1+g).*(1-X(:,1));
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,3);
            varargout = {f};
    end
end