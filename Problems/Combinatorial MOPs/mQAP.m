function varargout = mQAP(Operation,Global,input)
% <problem> <Combinatorial MOP>
% Multi-objective quadratic assignment problem
% c ---  0 --- Correlation parameter
% operator --- EApermutation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

persistent a b;

    c = Global.ParameterSet(0);
    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 10;
            Global.operator = @EApermutation;
            
            a = randi(100,Global.D,Global.D).*~eye(Global.D);
            r = arrayfun(@(s)rand(Global.D),1:Global.M,'UniformOutput',false);
            if c >= 0.9999
                [r{2:end}] = deal(r{1});
            elseif c >= 0
                r(2:end) = cellfun(@(S)norminv(S.*normcdf(1,r{1},1-sqrt(c))+(1-S).*normcdf(0,r{1},1-sqrt(c)),r{1},1-sqrt(c)),r(2:end),'UniformOutput',false);
            elseif c > -0.9999
                r(2:end) = cellfun(@(S)1-norminv(S.*normcdf(1,r{1},1-sqrt(-c))+(1-S).*normcdf(0,r{1},1-sqrt(-c)),r{1},1-sqrt(-c)),r(2:end),'UniformOutput',false);
            else
                [r{2:end}] = deal(1-r{1});
            end
            b     = cellfun(@(S)100*S.*~eye(Global.D),r,'UniformOutput',false);
            [a,b] = Global.ParameterFile(sprintf('mQAP-M%d-D%d-c%.4f',Global.M,Global.D,c),a,b);
            
            [~,PopDec] = sort(rand(input,Global.D),2);
            varargout  = {PopDec};
        case 'value'
            PopDec = input;
            M      = Global.M;
            
            PopObj = zeros(size(PopDec,1),M);
            [~,pi] = sort(PopDec,2);
            for i = 1 : M
                for j = 1 : size(PopDec,1)
                    PopObj(j,i) = sum(sum(a.*b{i}(pi(j,:),pi(j,:))));
                end
            end
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            RefPoint  = zeros(1,Global.M) + 10000*Global.D.^2;
            varargout = {RefPoint};
    end
end