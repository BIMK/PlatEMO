classdef MaF15 < PROBLEM
% <problem> <MaF>
% Inverted LSMOP8

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

    properties(Access = private)
        sublen;	% Number of variables in each subcomponent
        len;    % Cumulative sum of lengths of variable groups
    end
    methods
        %% Initialization
        function obj = MaF15()
            % Parameter setting
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 20*obj.Global.M;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = [ones(1,obj.Global.M-1),10.*ones(1,obj.Global.D-obj.Global.M+1)];
            obj.Global.encoding = 'real';
            % Calculate the number of variables in each subcomponent
            nk = 2;
            c  = 3.8*0.1*(1-0.1);
            for i = 1 : obj.Global.M-1
                c = [c,3.8.*c(end).*(1-c(end))];
            end
            obj.sublen = floor(c./sum(c).*obj.Global.D/nk);
            obj.len    = [0,cumsum(obj.sublen*nk)];
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D] = size(PopDec);
            M     = obj.Global.M;
            nk    = 2;
            PopDec(:,M:D) = (1+repmat(cos((M:D)./D*pi/2),N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
            for i = 1 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Griewank(PopDec(:,obj.len(i)+M-1+(j-1)*obj.sublen(i)+1:obj.len(i)+M-1+j*obj.sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : nk
                    G(:,i) = G(:,i) + Sphere(PopDec(:,obj.len(i)+M-1+(j-1)*obj.sublen(i)+1:obj.len(i)+M-1+j*obj.sublen(i)));
                end
            end
            G      = G./repmat(obj.sublen,N,1)./nk;
            PopObj = (1+G+[G(:,2:end),zeros(N,1)]).*(1-fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)]);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,obj.Global.M);
            P = 1 - P./repmat(sqrt(sum(P.^2,2)),1,obj.Global.M);
        end
    end
end

function f = Griewank(x)
    f = sum(x.^2,2)./4000 - prod(cos(x./repmat(sqrt(1:size(x,2)),size(x,1),1)),2) + 1;
end

function f = Sphere(x)
    f = sum(x.^2,2);
end