function varargout = ZDT5(Operation,Global,input)
% <problem> <ZDT>
% Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
% operator --- EAbinary

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
            Global.D        = 80;
            Global.D        = ceil(max(Global.D-30,1)/5)*5 + 30;
            Global.operator = @EAbinary;

            PopDec    = randi([0,1],input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            
            u      = zeros(size(PopDec,1),1+(size(PopDec,2)-30)/5);
            u(:,1) = sum(PopDec(:,1:30),2);
            for i = 2 : size(u,2)
                u(:,i) = sum(PopDec(:,(i-2)*5+31:(i-2)*5+35),2);
            end
            v           = zeros(size(u));
            v(u<5)      = 2 + u(u<5);
            v(u==5)     = 1;
            PopObj(:,1) = 1 + u(:,1);
            g           = sum(v(:,2:end),2);
            h           = 1./PopObj(:,1);
            PopObj(:,2) = g.*h;
            
            PopCon = [];
            
            varargout = {PopDec,PopObj,PopCon};
        case 'PF'
            f(:,1)    = 1 : 31;
            f(:,2)    = (Global.D-30)/5./f(:,1);
            varargout = {f};
    end
end