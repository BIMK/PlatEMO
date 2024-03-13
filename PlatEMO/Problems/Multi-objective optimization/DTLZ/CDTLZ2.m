classdef CDTLZ2 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% Convex DTLZ2

%------------------------------- Reference --------------------------------
% K. Deb and H. Jain, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part I:
% Solving problems with box constraints, IEEE Transactions on Evolutionary
% Computation, 2014, 18(4): 577-601.
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
            if isempty(obj.D); obj.D = obj.M + 9; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g      = sum((PopDec(:,obj.M:end)-0.5).^2,2);
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
            PopObj = [PopObj(:,1:obj.M-1).^4,PopObj(:,obj.M).^2];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R    = UniformPoint(N,obj.M).^2;
            temp = sum(sqrt(R(:,1:end-1)),2) + R(:,end);
            R    = R./[repmat(temp.^2,1,size(R,2)-1),temp];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                x = sin(a)*cos(a');
                y = sin(a)*sin(a');
                z = cos(a)*ones(size(a'));
                R = {x.^4,y.^4,z.^2};
            else
                R = [];
            end
        end
    end
end