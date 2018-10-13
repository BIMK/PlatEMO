function varargout = RMMEDA_F8(Operation,Global,input)
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
            Global.M        = 3;
            Global.M        = 3;
            Global.D        = 30;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            X = input;
            D = size(X,2);
            
            g = sum((X(:,3:end).^2-repmat(X(:,1),1,D-2)).^2,2);
            PopObj(:,1) = cos(pi/2*X(:,1)).*cos(pi/2*X(:,2)).*(1+g);
            PopObj(:,2) = cos(pi/2*X(:,1)).*sin(pi/2*X(:,2)).*(1+g);
            PopObj(:,3) = sin(pi/2*X(:,1)).*(1+g);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,3);
            f = f./repmat(sqrt(sum(f.^2,2)),1,3);
            varargout = {f};
    end
end