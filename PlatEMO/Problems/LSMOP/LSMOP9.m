classdef LSMOP9 < PROBLEM
% <problem> <LSMOP>
% Large-scale benchmark MOP
% nk --- 5 --- Number of subcomponents in each variable group

%------------------------------- Reference --------------------------------
% R. Cheng, Y. Jin, and M. Olhofer, Test problems for large-scale
% multiobjective and many-objective optimization, IEEE Transactions on
% Cybernetics, 2017, 47(12): 4108-4121.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        nk = 5; % Number of subcomponents in each variable group
        sublen;	% Number of variables in each subcomponent
        len;    % Cumulative sum of lengths of variable groups
    end
    methods
        %% Initialization
        function obj = LSMOP9()
            % Parameter setting
            obj.nk = obj.Global.ParameterSet(5);
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 100*obj.Global.M;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = [ones(1,obj.Global.M-1),10.*ones(1,obj.Global.D-obj.Global.M+1)];
            obj.Global.encoding = 'real';
            % Calculate the number of variables in each subcomponent
            c = 3.8*0.1*(1-0.1);
            for i = 1 : obj.Global.M-1
                c = [c,3.8.*c(end).*(1-c(end))];
            end
            obj.sublen = floor(c./sum(c).*obj.Global.D/obj.nk);
            obj.len    = [0,cumsum(obj.sublen*obj.nk)];
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D] = size(PopDec);
            M     = obj.Global.M;
            PopDec(:,M:D) = (1+repmat(cos((M:D)./D*pi/2),N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
            for i = 1 : 2 : M
                for j = 1 : obj.nk
                    G(:,i) = G(:,i) + Sphere(PopDec(:,obj.len(i)+M-1+(j-1)*obj.sublen(i)+1:obj.len(i)+M-1+j*obj.sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : obj.nk
                    G(:,i) = G(:,i) + Ackley(PopDec(:,obj.len(i)+M-1+(j-1)*obj.sublen(i)+1:obj.len(i)+M-1+j*obj.sublen(i)));
                end
            end
            G = 1 + sum(G./repmat(obj.sublen,N,1)./obj.nk,2);
            PopObj(:,1:M-1) = PopDec(:,1:M-1);
            PopObj(:,M)     = (1+G).*(M-sum(PopObj(:,1:M-1)./(1+repmat(G,1,M-1)).*(1+sin(3*pi.*PopObj(:,1:M-1))),2));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            interval     = [0,0.251412,0.631627,0.859401];
            median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
            X            = ReplicatePoint(N,obj.Global.M-1);
            X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
            X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
            P            = [X,2*(obj.Global.M-sum(X/2.*(1+sin(3*pi.*X)),2))];
        end
    end
end

function f = Sphere(x)
    f = sum(x.^2,2);
end

function f = Ackley(x)
    f = 20-20.*exp(-0.2.*sqrt(sum(x.^2,2)./size(x,2)))-exp(sum(cos(2.*pi.*x),2)./size(x,2))+exp(1);
end

function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end