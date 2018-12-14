classdef SDTLZ1 < PROBLEM
% <problem> <DTLZ variant>
% Scaled DTLZ1
% a --- 2 --- Scaling factor

%------------------------------- Reference --------------------------------
% K. Deb and H. Jain, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part I:
% Solving problems with box constraints, IEEE Transactions on Evolutionary
% Computation, 2014, 18(4): 577-601.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        a = 2;	% Scaling factor
    end
    methods
        %% Initialization
        function obj = SDTLZ1()
            obj.a = obj.Global.ParameterSet(2);
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = obj.Global.M + 4;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            M      = obj.Global.M;
            g      = 100*(D-M+1+sum((PopDec(:,M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,M:end)-0.5)),2));
            PopObj = 0.5*repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
            PopObj = PopObj.*repmat(obj.a.^(0:M-1),size(PopObj,1),1);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,obj.Global.M)/2;
            P = P.*repmat(obj.a.^(0:obj.Global.M-1),size(P,1),1);
        end
    end
end