classdef MaF13 < PROBLEM
% <problem> <MaF>
% P7

%------------------------------- Reference --------------------------------
% R. Cheng, M. Li, Y. Tian, X. Zhang, S. Yang, Y. Jin, and X. Yao, A
% benchmark test suite for evolutionary many-objective optimization,
% Complex & Intelligent Systems, 2017, 3(1): 67-81.
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
        function obj = MaF13()
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 5;
            end
            obj.Global.lower    = [zeros(1,2),zeros(1,obj.Global.D-2)-2];
            obj.Global.upper    = [ones(1,2),zeros(1,obj.Global.D-2)+2];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            [N,D] = size(X);
            Y = X - 2*repmat(X(:,2),1,D).*sin(2*pi*repmat(X(:,1),1,D)+repmat(1:D,N,1)*pi/D);
            PopObj(:,1) = sin(X(:,1)*pi/2)                   + 2*mean(Y(:,4:3:D).^2,2);
            PopObj(:,2) = cos(X(:,1)*pi/2).*sin(X(:,2)*pi/2) + 2*mean(Y(:,5:3:D).^2,2);
            PopObj(:,3) = cos(X(:,1)*pi/2).*cos(X(:,2)*pi/2) + 2*mean(Y(:,3:3:D).^2,2);
            PopObj(:,4:obj.Global.M) = repmat(PopObj(:,1).^2+PopObj(:,2).^10+PopObj(:,3).^10+2*mean(Y(:,4:D).^2,2),1,obj.Global.M-3);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,3);
            P = P./repmat(sqrt(sum(P.^2,2)),1,3);
            P = [P,repmat(P(:,1).^2+P(:,2).^10+P(:,3).^10,1,obj.Global.M-3)];
        end
    end
end