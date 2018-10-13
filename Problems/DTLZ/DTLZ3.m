function varargout = DTLZ3(Operation,Global,input)
% <problem> <DTLZ>
% Scalable Test Problems for Evolutionary Multi-Objective Optimization
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
            Global.D        = Global.M + 9;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,D]  = size(PopDec);
            M      = Global.M;
            
            g      = 100*(D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
            PopObj = repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,Global.M);
            f = f./repmat(sqrt(sum(f.^2,2)),1,Global.M);
            varargout = {f};
    end
end