classdef DTLZ6 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% Benchmark MOP proposed by Deb, Thiele, Laumanns, and Zitzler

%------------------------------- Reference --------------------------------
% K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, Scalable test problems
% for evolutionary multiobjective optimization, Evolutionary multiobjective
% Optimization. Theoretical Advances and Applications, 2005, 105-145.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+9; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g    = sum(PopDec(:,obj.M:end).^0.1,2);
            Temp = repmat(g,1,obj.M-2);
            PopDec(:,2:obj.M-1) = (1+2*Temp.*PopDec(:,2:obj.M-1))./(2+2*Temp);
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = [0:1/(N-1):1;1:-1/(N-1):0]';
            R = R./repmat(sqrt(sum(R.^2,2)),1,size(R,2));
            R = [R(:,ones(1,obj.M-2)),R];
            R = R./sqrt(2).^repmat([obj.M-2,obj.M-2:-1:0],size(R,1),1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M < 4
                R = obj.GetOptimum(100);
            else
                R = [];
            end
        end
    end
end