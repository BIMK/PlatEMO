classdef MaF13 < PROBLEM
% <multi/many> <real> <large/none>
% P7

%------------------------------- Reference --------------------------------
% R. Cheng, M. Li, Y. Tian, X. Zhang, S. Yang, Y. Jin, and X. Yao, A
% benchmark test suite for evolutionary many-objective optimization,
% Complex & Intelligent Systems, 2017, 3(1): 67-81.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
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
            if isempty(obj.D); obj.D = 5; end
            obj.M        = max(obj.M,3);
            obj.lower    = [zeros(1,2),zeros(1,obj.D-2)-2];
            obj.upper    = [ones(1,2),zeros(1,obj.D-2)+2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,D] = size(X);
            Y = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,N,1)*pi/D);
            PopObj(:,1) = sin(X(:,1)*pi/2)                   + 2*mean(Y(:,4:3:D).^2,2);
            PopObj(:,2) = cos(X(:,1)*pi/2).*sin(X(:,2)*pi/2) + 2*mean(Y(:,5:3:D).^2,2);
            PopObj(:,3) = cos(X(:,1)*pi/2).*cos(X(:,2)*pi/2) + 2*mean(Y(:,3:3:D).^2,2);
            PopObj(:,4:obj.M) = repmat(PopObj(:,1).^2+PopObj(:,2).^10+PopObj(:,3).^10+2*mean(Y(:,4:D).^2,2),1,obj.M-3);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,3);
            R = R./repmat(sqrt(sum(R.^2,2)),1,3);
            R = [R,repmat(R(:,1).^2+R(:,2).^10+R(:,3).^10,1,obj.M-3)];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
    end
end