classdef MinusIDTLZ3 < PROBLEM
% <2026> <multi/many> <real> <large/none> <expensive/none>
% Minus-Inverted DTLZ3

%------------------------------- Reference --------------------------------
% Y. Horaguchi, R. Kato, R. Nakahashi, J. Zhang, H. Ishibuchi, and M.
% Nakata. Evolutionary surrogate-assisted multiobjective optimization over
% high-dimensional search spaces: A comparative study. ACM Transactions on
% Evolutionary Learning and Optimization, 2026, 6(3): 25.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yuma Horaguchi and Ryudai Kato (email: horaguchi-yuma-hv@ynu.jp)

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
            g      = 100*(obj.D-obj.M+1+sum((PopDec(:,obj.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M:end)-0.5)),2));
            PopObj = (1+repmat(g,1,obj.M)) - repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(PopDec,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(PopDec,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
            PopObj = -PopObj;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = 1 - R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            g = (obj.D-obj.M+1)*220.2602;
            R = R*(1+g)*(-1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {1-sin(a)*cos(a'),1-sin(a)*sin(a'),1-cos(a)*ones(size(a'))};
                g = (obj.D-obj.M+1)*220.2602;
                for i = 1:obj.M
                    R{i} = R{i}*(1+g)*(-1);
                end
            else
                R = [];
            end
        end
    end
end