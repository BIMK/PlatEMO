classdef LSMOP4 < PROBLEM
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
        function obj = LSMOP4()
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
            PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
            G = zeros(N,M);
            for i = 1 : 2 : M
                for j = 1 : obj.nk
                    G(:,i) = G(:,i) + Ackley(PopDec(:,obj.len(i)+M-1+(j-1)*obj.sublen(i)+1:obj.len(i)+M-1+j*obj.sublen(i)));
                end
            end
            for i = 2 : 2 : M
                for j = 1 : obj.nk
                    G(:,i) = G(:,i) + Griewank(PopDec(:,obj.len(i)+M-1+(j-1)*obj.sublen(i)+1:obj.len(i)+M-1+j*obj.sublen(i)));
                end
            end
            G      = G./repmat(obj.sublen,N,1)./obj.nk;
            PopObj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,obj.Global.M);
        end
    end
end

function f = Ackley(x)
    f = 20-20.*exp(-0.2.*sqrt(sum(x.^2,2)./size(x,2)))-exp(sum(cos(2.*pi.*x),2)./size(x,2))+exp(1);
end

function f = Griewank(x)
    f = sum(x.^2,2)./4000 - prod(cos(x./repmat(sqrt(1:size(x,2)),size(x,1),1)),2) + 1;
end