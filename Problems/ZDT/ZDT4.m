function varargout = ZDT4(Operation,Global,input)
% <problem> <ZDT>
% Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
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
            Global.lower    = [0,zeros(1,Global.D-1)-5];
            Global.upper    = [1,zeros(1,Global.D-1)+5];
            Global.operator = @EAreal;
            

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            
            PopObj(:,1) = PopDec(:,1);
            g = 1+10*(size(PopDec,2)-1)+sum(PopDec(:,2:end).^2-10*cos(4*pi*PopDec(:,2:end)),2);
            h = 1-(PopObj(:,1)./g).^0.5;
            PopObj(:,2) = g.*h;
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^0.5;
            varargout = {f};
    end
end