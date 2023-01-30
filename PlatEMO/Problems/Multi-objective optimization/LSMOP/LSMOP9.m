classdef LSMOP9 < PROBLEM
% <multi/many> <real> <large/none>
% Large-scale benchmark MOP
% nk --- 5 --- Number of subcomponents in each variable group

%------------------------------- Reference --------------------------------
% R. Cheng, Y. Jin, and M. Olhofer, Test problems for large-scale
% multiobjective and many-objective optimization, IEEE Transactions on
% Cybernetics, 2017, 47(12): 4108-4121.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
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
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            obj.nk = obj.ParameterSet(5);
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = 100*obj.M; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = [ones(1,obj.M-1),10.*ones(1,obj.D-obj.M+1)];
            obj.encoding = ones(1,obj.D);
            % Calculate the number of variables in each subcomponent
            c = 3.8*0.1*(1-0.1);
            for i = 1 : obj.M-1
                c = [c,3.8.*c(end).*(1-c(end))];
            end
            obj.sublen = floor(c./sum(c).*(obj.D-obj.M+1)/obj.nk);
            obj.len    = [0,cumsum(obj.sublen*obj.nk)];
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D] = size(PopDec);
            M     = obj.M;
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
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            interval     = [0,0.251412,0.631627,0.859401];
            median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
            X            = UniformPoint(N,obj.M-1,'grid');
            X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
            X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
            R            = [X,2*(obj.M-sum(X/2.*(1+sin(3*pi.*X)),2))];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                x      = linspace(0,1,100)';
                y      = 2*(2-x/2.*(1+sin(3*pi*x)));
                nd     = NDSort([x,y],1)==1;
                x(~nd) = nan;
                R      = [x,y];
            elseif obj.M == 3
                [x,y]  = meshgrid(linspace(0,1,20));
                z      = 2*(3-x/2.*(1+sin(3*pi*x))-y/2.*(1+sin(3*pi*y)));
                nd     = reshape(NDSort([x(:),y(:),z(:)],1)==1,size(z));
                z(~nd) = nan;
                R      = {x,y,z};
            else
                R = [];
            end
        end
    end
end

function f = Sphere(x)
    f = sum(x.^2,2);
end

function f = Ackley(x)
    f = 20-20.*exp(-0.2.*sqrt(sum(x.^2,2)./size(x,2)))-exp(sum(cos(2.*pi.*x),2)./size(x,2))+exp(1);
end