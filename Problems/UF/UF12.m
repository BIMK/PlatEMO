function varargout = UF12(Operation,Global,input)
% <problem> <UF>
% Multiobjective optimization Test Instances for the CEC 2009 Special
% Session and Competition
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

persistent Bound Lamda M;

    % This problem is R2_DTLZ3_M5, the numbers of objectives and decision
    % variables are fixed to 5 and 30, respectively
    switch Operation
        case 'init'
            load('UF12_parameter','Bound','Lamda','M');
            
            Global.M        = 5;
            Global.M        = 5;
            Global.D        = 30;
            Global.D        = 30;
            Global.lower    = Bound(1,:);
            Global.upper    = Bound(2,:);
            Global.operator = @EAreal;

            PopDec    = rand(input,size(Bound,2)).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            
            lamda = repmat(Lamda,size(PopDec,1),1);
            z     = PopDec*M';
            p     = zeros(size(z));
            temp1 = z < 0;
            temp2 = z > 1;
            p(temp1) = -z(temp1);
            p(temp2) = z(temp2) - 1;
            z(temp1) = -lamda(temp1).*z(temp1);
            z(temp2) = 1 - lamda(temp2).*(z(temp2)-1);
            psum = zeros(size(PopDec,1),5);
            for i = 1 : 5
                psum(:,i) = sqrt(sum(p(:,[1:min(6-i,4),5:end]).^2,2));
            end
            g      = 100*(26+sum((z(:,5:end)-0.5).^2-cos(20.*pi.*(z(:,5:end)-0.5)),2));
            PopObj = 2./(1+exp(-psum)).*(1+repmat(1+g,1,5).*fliplr(cumprod([ones(size(g,1),1),cos(z(:,1:5-1)*pi/2)],2)).*[ones(size(g,1),1),sin(z(:,5-1:-1:1)*pi/2)]);

            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = UniformPoint(input,5);
            f = f./repmat(sqrt(sum(f.^2,2)),1,5) + 1;
            varargout = {f};
    end
end