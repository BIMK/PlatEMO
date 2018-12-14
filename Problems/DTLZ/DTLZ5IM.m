classdef DTLZ5IM < PROBLEM
% <problem> <DTLZ variant>
% Extended DTLZ5
% I --- 2 --- Dimensionality of the true Pareto front

%------------------------------- Reference --------------------------------
% K. Deb and D. K. Saxena, On finding Pareto-optimal solutions through
% dimensionality reduction for certain large-dimensional multi-objective
% optimization problems, KanGAL Report 2005011, 2005.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        I = 2;	% Dimensionality of the true Pareto front
    end
    methods
        %% Initialization
        function obj = DTLZ5IM()
            obj.I = obj.Global.ParameterSet(2);
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
            M      = obj.Global.M;
            g      = sum((PopDec(:,M:end)-0.5).^2,2);
            Temp   = repmat(g,1,M-obj.I);
            PopDec(:,obj.I:M-1) = (1+2*Temp.*PopDec(:,obj.I:M-1))./(2+2*Temp);
            PopObj = repmat(1+100*g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)];
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,obj.I);
            P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
            P = [P(:,ones(1,obj.Global.M-size(P,2))),P];
            P = P./sqrt(2).^repmat(max([obj.Global.M-obj.I,obj.Global.M-obj.I:-1:2-obj.I],0),size(P,1),1);
        end
    end
end