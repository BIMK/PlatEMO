classdef SMD10 < PROBLEM
% <multi> <real> <constrained> <bilevel>
% Bilevel optimization problems proposed by Sinha, Malo, and Deb
% maxFElower --- 500 --- Maximum number of lower level function evaluations for each solution

%------------------------------- Reference --------------------------------
% A. Sinha, P. Malo, K. Deb, Test problem construction for single-objective 
% bilevel optimization, Evolutionary Computation, 2014, 22(3): 439-477.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    properties(SetAccess = private)
        maxFElower; % Maximum number of lower level function evaluations for each solution
        DU;         % Number of decision variables of the upper level
        DL;         % Number of decision variables of the lower level
        C;       	% Number of upper constraints
        p;          % The length of xu1
        q;          % The length of xl1
        r;          % The length of xu2 and xl2
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.maxFElower = obj.ParameterSet(500);
            obj.M = 2;
            if isempty(obj.D); obj.D = 5; end
            obj.DU = floor(obj.D/2);
            obj.DL = obj.D - obj.DU;
            obj.C  = obj.DU;
            obj.r  = floor(obj.DU/2);
            obj.p  = obj.DU - obj.r;
            obj.q  = obj.DL - obj.r;
            obj.lower    = [-5*ones(1,obj.DU), -5*ones(1,obj.q), -pi/2*ones(1,obj.r)+1e-6];
            obj.upper    = [10*ones(1,obj.DU), 10*ones(1,obj.q),  pi/2*ones(1,obj.r)-1e-6];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate upper level and lower level objective values
        function PopObj = CalObj(obj,PopDec)
            xu1 = PopDec(:,1:obj.p); 
            xu2 = PopDec(:,obj.p+1:obj.p+obj.r); 
            xl1 = PopDec(:,obj.p+obj.r+1:obj.p+obj.r+obj.q); 
            xl2 = PopDec(:,obj.p+obj.r+obj.q+1:end);
            % Calculate  UL objective values
            PopObj(:,1) = sum((xu1-2).^2,2) + sum(xl1.^2,2) + sum((xu2-2).^2,2) - sum((xu2-tan(xl2)).^2,2);
            % Calculate  LL objective values
            PopObj(:,2) = sum(xu1.^2,2) + sum((xl1-2).^2,2) + sum((xu2-tan(xl2)).^2,2);
        end
        %% Calculate upper level and lower level constraint violations
        function PopCon = CalCon(obj,PopDec)
            xu1 = PopDec(:,1:obj.p); 
            xu2 = PopDec(:,obj.p+1:obj.p+obj.r); 
            xl1 = PopDec(:,obj.p+obj.r+1:obj.p+obj.r+obj.q); 
            xl2 = PopDec(:,obj.p+obj.r+obj.q+1:end);
            % Upper level constraint violation
            PopCon(:,1:obj.p) = -xu1 - xu1.^3 + sum(xu1.^3,2) + sum(xu2.^3,2);
            PopCon(:,obj.p+1:obj.p+obj.r) = -xu2 - xu2.^3 + sum(xu2.^3,2) + sum(xu1.^3,2);
            % Lower level constraint violation
            PopCon(:,obj.p+obj.r+1:obj.p+obj.r+obj.q) = -xl1 - xl1.^3 + sum(xl1.^3,2);
        end
        %% Calculate lower level objective values
        function llPopulation = EvaluationLower(obj,varargin)
            PopDec            = obj.CalDec(varargin{1});
            PopObj            = obj.CalObj(PopDec);
            PopObj(:,1)       = nan;
            PopCon            = obj.CalCon(PopDec);
            PopCon(:,1:obj.C) = nan;
            llPopulation      = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
        end
    end
end