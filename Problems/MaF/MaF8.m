function varargout = MaF8(Operation,Global,input)
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

persistent Points;

    % This problem is multi-point distance minimization problem (similar to
    % Pareto-box problem but more difficult)
    switch Operation
        case 'init'
            Global.M          = 10;
            Global.D          = 2;
            Global.D          = 2;
            Global.lower      = [-10000,-10000];
            Global.upper      = [10000,10000];
            Global.operator   = @EAreal;
            Global.evaluation = max(1e5,1e4*Global.D);
            
            Points = [];
            [thera,rho] = cart2pol(0,1);
            [Points(:,1),Points(:,2)] = pol2cart(thera-(1:Global.M)*2*pi/Global.M,rho);
            
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            
            PopObj = pdist2(PopDec,Points);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            [X,Y]     = ndgrid(linspace(-1,1,ceil(sqrt(input))));
            ND        = inpolygon(X(:),Y(:),Points(:,1),Points(:,2));
            PopObj    = pdist2([X(ND),Y(ND)],Points);
            varargout = {PopObj};
        case 'draw'
            cla; Draw(input);
            plot(Points([1:end,1],1),Points([1:end,1],2),'-k','LineWidth',1.5);
            plot(Points(:,1),Points(:,2),'ok','MarkerSize',6,'Marker','o','Markerfacecolor',[1 1 1],'Markeredgecolor',[.4 .4 .4]);
            xlabel('\itx\rm_1'); ylabel('\itx\rm_2');
    end
end