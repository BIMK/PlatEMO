classdef MaF2 < PROBLEM
% <multi/many> <real> <large/none>
% DTLZ2BZ

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
            if isempty(obj.D); obj.D = obj.M + 9; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g = zeros(size(PopDec,1),obj.M);
            for m = 1 : obj.M
                if m < obj.M
                    g(:,m) = sum(((PopDec(:,obj.M+(m-1)*floor((obj.D-obj.M+1)/obj.M):obj.M+m*floor((obj.D-obj.M+1)/obj.M)-1)/2+1/4)-0.5).^2,2);
                else
                    g(:,m) = sum(((PopDec(:,obj.M+(obj.M-1)*floor((obj.D-obj.M+1)/obj.M):obj.D)/2+1/4)-0.5).^2,2);
                end
            end
            PopObj = (1+g).*fliplr(cumprod([ones(size(g,1),1),cos((PopDec(:,1:obj.M-1)/2+1/4)*pi/2)],2)).*[ones(size(g,1),1),sin((PopDec(:,obj.M-1:-1:1)/2+1/4)*pi/2)];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            c = zeros(size(R,1),obj.M-1);
            for i = 1 : size(R,1)
                for j = 2 : obj.M
                    temp = R(i,j)/R(i,1)*prod(c(i,obj.M-j+2:obj.M-1));
                    c(i,obj.M-j+1) = sqrt(1/(1+temp^2));
                end
            end
            if obj.M > 5
                c = c.*(cos(pi/8)-cos(3*pi/8)) + cos(3*pi/8);
            else
                c(any(c<cos(3*pi/8)|c>cos(pi/8),2),:) = [];
            end
            R = fliplr(cumprod([ones(size(c,1),1),c(:,1:obj.M-1)],2)).*[ones(size(c,1),1),sqrt(1-c(:,obj.M-1:-1:1).^2)];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                x      = linspace(0,pi/2,100)';
                R(:,1) = cos(x);
                R(:,2) = sin(x);
                R(cos(x)<cos(3*pi/8)|cos(x)>cos(pi/8),:) = nan;
            elseif obj.M == 3
                a  = linspace(0,pi/2,20)';
                x  = cos(a)*cos(a');
                y  = cos(a)*sin(a');
                z  = sin(a)*ones(size(a'));
                c1 = cos(a*ones(size(a')));
                c2 = cos(ones(size(a))*a');
                z(c1<cos(3*pi/8)|c2<cos(3*pi/8)|c1>cos(pi/8)|c2>cos(pi/8)) = nan;
                R = {x,y,z};
            else
                R = [];
            end
        end
    end
end