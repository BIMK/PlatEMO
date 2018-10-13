function varargout = MONRP(Operation,Global,input)
% <problem> <Combinatorial MOP>
% Multi-objective Next Release Problem
% m --- 100 --- Number of customers
% operator  --- EAbinary

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

persistent cost value;

    m = Global.ParameterSet(100);
    switch Operation
        case 'init'
            Global.M        = 2;
            Global.M        = 2;
            Global.D        = 100;
            Global.operator = @EAbinary;
            
            n            = Global.D;
            [cost,value] = Global.ParameterFile(sprintf('MONRP-n%d-m%d',n,m),randi(9,1,n),randi([0 5],n,m));
            
            PopDec    = randi([0,1],input,n);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            
            PopObj(:,1) = sum(repmat(cost,size(PopDec,1),1).*PopDec,2);
            PopObj(:,2) = sum(value(:)) - sum(PopDec*value,2);
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            RefPoint  = [sum(cost),sum(value(:))];
            varargout = {RefPoint};
        case 'draw'
            PopObj(:,1) = sum(repmat(cost,size(input,1),1).*input,2);
            PopObj(:,2) = sum(input*value,2);
            cla; Draw(PopObj);
            xlabel('Cost'); ylabel('Satisfaction Score');
    end
end