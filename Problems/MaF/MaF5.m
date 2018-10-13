function varargout = MaF5(Operation,Global,input)
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

    % This problem is scaled DTLZ4
    switch Operation
        case 'init'
            Global.M          = 3;
            Global.D          = Global.M + 9;
            Global.lower      = zeros(1,Global.D);
            Global.upper      = ones(1,Global.D);
            Global.operator   = @EAreal;
            Global.evaluation = max(1e5,1e4*Global.D);

            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            M      = Global.M;

            PopDec(:,1:M-1) = PopDec(:,1:M-1).^100;
            g      = sum((PopDec(:,M:end)-0.5).^2,2);
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
            PopObj = PopObj.*repmat(2.^(M:-1:1),size(g,1),1);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,Global.M);
            f = f./repmat(sqrt(sum(f.^2,2)),1,Global.M);
            f = f.*repmat(2.^(Global.M:-1:1),size(f,1),1);
            varargout = {f};
    end
end