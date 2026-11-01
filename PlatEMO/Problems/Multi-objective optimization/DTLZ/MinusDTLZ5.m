classdef MinusDTLZ5 < PROBLEM
% <2017> <multi/many> <real> <large/none> <expensive/none>
% Minus-DTLZ5

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
            g    = sum((PopDec(:,obj.M:end)-0.5).^2,2);
            Temp = repmat(g,1,obj.M-2);
            PopDec(:,2:obj.M-1) = (1+2*Temp.*PopDec(:,2:obj.M-1))./(2+2*Temp);
            PopObj = -1*repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            g = 0.5^2*(obj.D-obj.M+1);
            cosmin = cos(pi / (4 * (1 + g)));
            cosmax = cos(pi / (4 * (1 + g)) * (1 + 2 * g));

            c = zeros(size(R,1),obj.M-1);
            for i = 1 : size(R,1)
                for j = 2 : obj.M
                    temp = R(i,j)/R(i,1)*prod(c(i,obj.M-j+2:obj.M-1)); % tan
                    c(i,obj.M-j+1) = sqrt(1/(1+temp^2));               % cos
                end
            end
            if obj.M > 5
                c(:, 2 : end) = c(:, 2 : end).*(cosmin-cosmax) + cosmax;
            else
                c(any(c(:, 2 : end)<cosmax|c(:, 2 : end)>cosmin,2),:) = [];
            end
            R = fliplr(cumprod([ones(size(c,1),1),c(:,1:obj.M-1)],2)).*[ones(size(c,1),1),sqrt(1-c(:,obj.M-1:-1:1).^2)];
            R = R*(1+g)*(-1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                g = 0.5^2*(obj.D-obj.M+1);
                a1 = linspace(0,pi/2,20)';
                a2 = linspace(pi / (4 * (1 + g)),pi / (4 * (1 + g)) * (1 + 2 * g),20)';
                x  = cos(a1)*cos(a2');
                y  = cos(a1)*sin(a2');
                z  = sin(a1)*ones(size(a1'));
                R = {-(1+g)*x,-(1+g)*y,-(1+g)*z};
            else
                R = [];
            end
        end
    end
end