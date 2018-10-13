function varargout = MaF6(Operation,Global,input)
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

    % This problem is DTLZ5(I,M) with I=2
    I = 2;
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
            Temp   = repmat(g,1,M-I);
            PopDec(:,I:M-1) = (1+2*Temp.*PopDec(:,I:M-1))./(2+2*Temp);
            PopObj = repmat(1+100*g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,I);
            f = f./repmat(sqrt(sum(f.^2,2)),1,size(f,2));
            f = [f(:,ones(1,Global.M-size(f,2))),f];
            f = f./sqrt(2).^repmat(max([Global.M-I,Global.M-I:-1:2-I],0),size(f,1),1);
            varargout = {f};
    end
end