function varargout = MaF1(Operation,Global,input)
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

    % This problem is inverted DTLZ1 with the g function of DTLZ2
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

            g      = sum((PopDec(:,M:end)-0.5).^2,2);
            PopObj = repmat(1+g,1,M) - repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),PopDec(:,1:M-1)],2)).*[ones(size(g,1),1),1-PopDec(:,M-1:-1:1)];
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,Global.M);
            varargout = {1-f};
    end
end