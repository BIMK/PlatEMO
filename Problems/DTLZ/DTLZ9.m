function varargout = DTLZ9(Operation,Global,input)
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
            Global.M        = 2;
            Global.D        = 10*Global.M;
            Global.D        = ceil(Global.D/Global.M)*Global.M;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input.^0.1;
            [N,D]  = size(PopDec);
            M      = Global.M;
            
            PopObj = zeros(N,M);
            for m = 1 : M
                PopObj(:,m) = sum(PopDec(:,(m-1)*D/M+1:m*D/M),2);
            end
            
            PopCon = 1 - repmat(PopObj(:,M).^2,1,M-1) - PopObj(:,1:M-1).^2;
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            Temp      = (0:1/(input-1):1)';
            f         = [repmat(cos(0.5.*pi.*Temp),1,Global.M-1),sin(0.5.*pi.*Temp)];
            varargout = {f};
    end
end