function varargout = MOEADM2M_F5(Operation,Global,input)
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
            Global.M        = 2;
            Global.M        = 2;
            Global.D        = 10;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            X = input;
            
            t = X(:,2:end) - repmat(sin(pi/2*X(:,1)),1,size(X,2)-1);
            g = 2*abs(cos(pi*X(:,1))).*sum(-0.9*t.^2+abs(t).^0.6,2);
            PopObj(:,1) = (1+g).*X(:,1);
            PopObj(:,2) = (1+g).*(1-sqrt(X(:,1)));
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1 - sqrt(f(:,1));
            varargout = {f};
    end
end