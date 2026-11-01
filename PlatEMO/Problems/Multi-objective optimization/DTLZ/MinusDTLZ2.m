classdef MinusDTLZ2 < PROBLEM
% <2017> <multi/many> <real> <large/none> <expensive/none>
% Minus-DTLZ2

%------------------------------- Reference --------------------------------
% H. Ishibuchi, Y. Setoguchi, H. Masuda, and Y. Nojima. Performance of
% decomposition-based many-objective algorithms strongly depends on Pareto
% front shapes. IEEE Transactions on Evolutionary Computation, 2017, 21(2):
% 169-190.
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
            g      = sum((PopDec(:,obj.M:end)-0.5).^2,2);
            PopObj = -repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            g = 0.5^2*(obj.D-obj.M+1);
            R = R*(1+g)*(-1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
                g = 0.5^2*(obj.D-obj.M+1);
                for i = 1:obj.M
                    R{i} = R{i}*(1+g)*(-1);
                end
            else
                R = [];
            end
        end
    end
end