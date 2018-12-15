classdef MaF2 < PROBLEM
% <problem> <MaF>
% DTLZ2BZ

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
        function obj = MaF2()
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = obj.Global.M + 9;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D] = size(PopDec);
            M     = obj.Global.M;
            g     = zeros(N,M);
            for m = 1 : M
                if m < M
                    g(:,m) = sum(((PopDec(:,M+(m-1)*floor((D-M+1)/M):M+m*floor((D-M+1)/M)-1)/2+1/4)-0.5).^2,2);
                else
                    g(:,m) = sum(((PopDec(:,M+(M-1)*floor((D-M+1)/M):D)/2+1/4)-0.5).^2,2);
                end
            end
            PopObj = (1+g).*fliplr(cumprod([ones(size(g,1),1),cos((PopDec(:,1:M-1)/2+1/4)*pi/2)],2)).*[ones(size(g,1),1),sin((PopDec(:,M-1:-1:1)/2+1/4)*pi/2)];
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            M = obj.Global.M;
            P = UniformPoint(N,M);
            c = zeros(size(P,1),M-1);
            for i = 1 : size(P,1)
                for j = 2 : M
                    temp = P(i,j)/P(i,1)*prod(c(i,M-j+2:M-1));
                    c(i,M-j+1) = sqrt(1/(1+temp^2));
                end
            end
            if M > 5
                c = c.*(cos(pi/8)-cos(3*pi/8)) + cos(3*pi/8);
            else
                c(any(c<cos(3*pi/8)|c>cos(pi/8),2),:) = [];
            end
            P = fliplr(cumprod([ones(size(c,1),1),c(:,1:M-1)],2)).*[ones(size(c,1),1),sqrt(1-c(:,M-1:-1:1).^2)];
        end
    end
end