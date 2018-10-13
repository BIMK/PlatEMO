function varargout = RMMEDA_F1(Operation,Global,input)
% <problem> <RMMEDA>
% RM-MEDA: A Regularity Model-Based Multiobjective Estimation of
% Distribution Algorithm
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
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            X = input;
            D = size(X,2);
            
            g = 1 + 9*mean((X(:,2:end)-repmat(X(:,1),1,D-1)).^2,2);
            PopObj(:,1) = X(:,1);
            PopObj(:,2) = g.*(1-sqrt(PopObj(:,1)./g));
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1 - sqrt(f(:,1));
            varargout = {f};
    end
end