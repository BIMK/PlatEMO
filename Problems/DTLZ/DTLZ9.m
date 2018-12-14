classdef DTLZ9 < PROBLEM
% <problem> <DTLZ>
% Benchmark MOP proposed by Deb, Thiele, Laumanns, and Zitzler

%------------------------------- Reference --------------------------------
% K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, Scalable test problems
% for evolutionary multiobjective optimization, Evolutionary multiobjective
% Optimization. Theoretical Advances and Applications, 2005, 105-145.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function obj = DTLZ9()
            if isempty(obj.Global.M)
                obj.Global.M = 2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 10*obj.Global.M;
            end
            obj.Global.D        = ceil(obj.Global.D/obj.Global.M)*obj.Global.M;
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            PopDec = PopDec.^0.1;
            M      = obj.Global.M;
            PopObj = zeros(N,M);
            for m = 1 : M
                PopObj(:,m) = sum(PopDec(:,(m-1)*D/M+1:m*D/M),2);
            end
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            M      = obj.Global.M;
            PopObj = obj.CalObj(PopDec);
            PopCon = 1 - repmat(PopObj(:,M).^2,1,M-1) - PopObj(:,1:M-1).^2;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            Temp = (0:1/(N-1):1)';
            P    = [repmat(cos(0.5.*pi.*Temp),1,obj.Global.M-1),sin(0.5.*pi.*Temp)];
        end
    end
end