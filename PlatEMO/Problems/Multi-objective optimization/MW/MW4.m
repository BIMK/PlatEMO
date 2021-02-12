classdef MW4 < PROBLEM
% <multi/many> <real> <large/none> <constrained>
% Constrained benchmark MOP proposed by Ma and Wang

%------------------------------- Reference --------------------------------
% Z. Ma and Y. Wang, Evolutionary constrained multiobjective optimization:
% Test suite construction and performance comparisons. IEEE Transactions on
% Evolutionary Computation, 2019, 23(6): 972-986.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
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
            if isempty(obj.D); obj.D = 15; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
             g      = sum(1 - exp(-10*((X(:,obj.M:end).^(obj.D-obj.M)) - 0.5 - ((obj.M:obj.D) - 1)/(2*obj.D)).^2),2);
             PopObj = repmat(1+g,1,obj.M).*flip(cumprod([ones(size(X,1),1),X(:,1:obj.M-1)],2),2).*[ones(size(X,1),1),1-X(:,obj.M-1:-1:1)];
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj = obj.CalObj(X);
            l      = PopObj(:,end) - sum(PopObj(:,1:(end-1)),2);
            PopCon = sum(PopObj,2) - (1 + 0.4*sin(2.5*pi*l).^8);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            l = R(:,end) - sum(R(:,1:end-1),2);
            c = (1+0.4*sin(2.5*pi*l).^8) - sum(R,2);
            R(c<0,:) = [];
        end
        %% Generate the feasible region
        function R = GetPF(obj)
            if obj.M == 2
                [x,y] = meshgrid(linspace(0,1.2,400));
                z     = nan(size(x));
                fes   = x + y - 1 - 0.4*sin(2.5*pi*(y-x)).^8 <= 0;
                z(fes & x+y>=1) = 0;
                R = {x,y,z};
            elseif obj.M == 3
                a = linspace(0,1,10)';
                x = a*a';
                y = a*(1-a');
                z = (1-a)*ones(size(a'));
                R = {x,y,z};
            else
                R = [];
            end
        end
    end
end