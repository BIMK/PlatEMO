classdef PROBLEM < handle & matlab.mixin.Heterogeneous
%PROBLEM - The superclass of problems.
%
%   This is the superclass of problems. An object of PROBLEM stores all the
%   settings of the problem.
%
% PROBLEM properties:
%   N               <read-only> population size
%   M               <read-only> number of objectives
%   D               <read-only> number of decision variables
%   maxFE           <read-only> maximum number of function evaluations
%   FE              <read-only> number of consumed function evaluations
%   encoding        <read-only> encoding scheme of decision variables
%   lower           <read-only> lower bound of decision variables
%   upper           <read-only> upper bound of decision variables
%   optimum         <read-only> optimal values of the problem
%   PF              <read-only> image of Pareto front
%   parameter       <read-only> other parameters of the problem
%
% PROBLEM methods:
%   PROBLEM         <protected> the constructor, which sets all the properties specified by user
%   Setting         <public>    default settings of the problem
%   Initialization 	<public>    generate initial solutions
%   CalDec          <public>    repair invalid solutions
%   CalObj          <public>    calculate the objective values of solutions
%   CalCon          <public>    calculate the constraint violations of solutions
%   GetOptimum      <public>    generate the optimums of the problem
%   GetPF          	<public>    generate the image of Pareto front
%   DrawDec         <public>    display a population in the decision space
%   DrawObj         <public>    display a population in the objective space
%   Current         <static>    get or set the current PROBLEM object
%   ParameterSet	<protected>	obtain the parameters of the problem

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        N        = 100;             	% Population size
        FE       = 0;                   % Number of consumed function evaluations
    end
    properties(SetAccess = protected)
        M;                             	% Number of objectives
        D;                             	% Number of decision variables
        maxFE    = 10000;              	% Maximum number of function evaluations
        encoding = 'real';            	% Encoding scheme of decision variables
        lower    = 0;                 	% Lower bound of decision variables
        upper    = 1;                  	% Upper bound of decision variables
        optimum;                       	% Optimal values of the problem
        PF;                          	% Image of Pareto front
        parameter = {};                	% Other parameters of the problem
    end
    methods(Access = protected)
        function obj = PROBLEM(varargin)
        %PROBLEM - The constructor of PROBLEM.
        %
        %   Problem = proName('Name',Value,'Name',Value,...) generates an
        %   object with the properties specified by the inputs. proName is
        %   PROBLEM or a subclass of PROBLEM.
        %
        %   If proName is PROBLEM, the properties objFcn and conFcn should
        %   be specified to define the objective and constraint functions.
        %   If proName is a subclass of PROBLEM, the properties objFcn and
        %   conFcn are useless since the objective and constraint functions
        %   are defined in the methods proName.CalObj and proName.CalCon.
        %
        %   The properties M, D, maxFE, encoding, lower, and upper will be
        %   further revised in the method proName.Setting.
        %
        %   Example:
        %       Problem = PROBLEM('objFcn',@(x)sum(x,2))
        %       Problem = DTLZ2('M',5,'D',10)

            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'N','M','D','maxFE','parameter'}))
                obj.(varargin{i}) = varargin{i+1};
            end
            obj.Setting();
            obj.optimum = obj.GetOptimum(10000);
            obj.PF      = obj.GetPF();
        end
    end
    methods
        function Setting(obj)
        %Setting - Default settings of the problem.
        %
        %   This function is expected to be implemented in each subclass of
        %   PROBLEM, which will be called automatically.
        end
        function Population = Initialization(obj,N)
        %Initialization - Generate initial solutions.
        %
        %   P = obj.Initialization() randomly generates the decision
        %   variables of obj.N solutions and returns the SOLUTION objects.
        %
        %   P = obj.Initialization(N) generates N solutions.
        %
        %   Example:
        %       Population = Problem.Initialization()
        
            if nargin < 2
            	N = obj.N;
            end
            switch obj.encoding
                case 'binary'
                    PopDec = randi([0,1],N,obj.D);
                case 'permutation'
                    [~,PopDec] = sort(rand(N,obj.D),2);
                otherwise
                    PopDec = unifrnd(repmat(obj.lower,N,1),repmat(obj.upper,N,1));
            end
            Population = SOLUTION(PopDec);
        end
        function PopDec = CalDec(obj,PopDec)
        %CalDec - Repair invalid solutions.
        %
        %   Dec = obj.CalDec(Dec) repairs the invalid (not infeasible)
        %   decision variables of Dec.
        %
        %   An invalid solution indicates that it is out of the decision
        %   space, while an infeasible solution indicates that it does not
        %   satisfy all the constraints.
        %
        %   This function will be used when obj.decFcn is not specified.
        %
        %   Example:
        %       PopDec = Problem.CalDec(PopDec)

            if strcmp(obj.encoding,'real')
                PopDec = max(min(PopDec,repmat(obj.upper,size(PopDec,1),1)),repmat(obj.lower,size(PopDec,1),1));
            end
        end
        function PopObj = CalObj(obj,PopDec)
        %CalObj - Calculate the objective values of solutions.
        %
        %   Obj = obj.CalObj(Dec) returns the objective values of Dec.
        %
        %   This function will be used when obj.objFcn is not specified.
        %
        %   Example:
        %       PopObj = Problem.CalObj(PopDec)

            PopObj = zeros(size(PopDec,1),obj.M);
        end
        function PopCon = CalCon(obj,PopDec)
        %CalCon - Calculate the constraint violations of solutions.
        %
        %   Con = obj.CalCon(Dec) returns the constraint violations of Dec.
        %
        %   This function will be used when obj.conFcn is not specified.
        %
        %   Example:
        %       PopCon = Problem.CalCon(PopDec)
        
            PopCon = zeros(size(PopDec,1),1);
        end
        function R = GetOptimum(obj,N)
        %GetOptimum - Generate the optimums of the problem.
        %
        %   R = obj.GetOptimum(N) returns N optimums of the problem for
        %   metric calculation.
        %
        %   For multi-objective optimization problems, an optimum can be a
        %   point on the Pareto front; if the Pareto front is unknown, an
        %   optimum can be a reference point for hypervolume calculation.
        %
        %   For multi-modal multi-objective optimization problems, an
        %   optimum can be the decision variables of a Pareto optimal
        %   solution.
        %
        %   For single-objective optimization problems, an optimum can be
        %   the minimum objective value of the problem.
        %
        %   Example:
        %       R = Problem.GetOptimum(10000)
        
            if obj.M > 1
                R = ones(1,obj.M);
            else
                R = 0;
            end
        end
        function R = GetPF(obj)
        %GetPF - Generate the image of Pareto front.
        %
        %   R = obj.GetPF() returns the image of Pareto front for objective
        %   visualization.
        %
        %   For bi-objective optimization problems, the image should be a
        %   one-dimensional curve.
        %
        %   For tri-objective optimization problems, the image should be a
        %   two-dimensional surface.
        %
        %   For constrained multi-objective optimization problems, the
        %   image can be the feasible region.
        %
        %   Example:
        %       R = Problem.GetPF()
        
            R = [];
        end
        function DrawDec(obj,Population)
        %DrawDec - Display a population in the decision space.
        %
        %   obj.DrawDec(P) displays the decision variables of population P.
        %
        %   Example:
        %       Problem.DrawDec(Population)
        
            if strcmp(obj.encoding,'binary')
                Draw(logical(Population.decs));
            else
                Draw(Population.decs,{'\it x\rm_1','\it x\rm_2','\it x\rm_3'});
            end
        end
        function DrawObj(obj,Population)
        %DrawObj - Display a population in the objective space.
        %
        %   obj.DrawObj(P) displays the objective values of population P.
        %
        %   Example:
        %       Problem.DrawObj(Population)

            ax = Draw(Population.objs,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
            if ~isempty(obj.PF)
                if ~iscell(obj.PF)
                    if obj.M == 2
                        plot(ax,obj.PF(:,1),obj.PF(:,2),'-k','LineWidth',1);
                    elseif obj.M == 3
                        plot3(ax,obj.PF(:,1),obj.PF(:,2),obj.PF(:,3),'-k','LineWidth',1);
                    end
                else
                    if obj.M == 2
                        surf(ax,obj.PF{1},obj.PF{2},obj.PF{3},'EdgeColor','none','FaceColor',[.85 .85 .85]);
                    elseif obj.M == 3
                        surf(ax,obj.PF{1},obj.PF{2},obj.PF{3},'EdgeColor',[.8 .8 .8],'FaceColor','none');
                    end
                    set(ax,'Children',ax.Children(flip(1:end)));
                end
            elseif size(obj.optimum,1) > 1 && obj.M < 4
                if obj.M == 2
                    plot(ax,obj.optimum(:,1),obj.optimum(:,2),'.k');
                elseif obj.M == 3
                    plot3(ax,obj.optimum(:,1),obj.optimum(:,2),obj.optimum(:,3),'.k');
                end
            end
        end
    end
    methods(Static, Sealed)
        function obj = Current(obj)
        %Current - Get or set the current PROBLEM object.
        %
        %   Pro = PROBLEM.Current() returns the current PROBLEM object.
        %
        %   PROBLEM.Current(Pro) sets the current PROBLEM object to Pro and
        %   sets Pro.evaluated to 0.
        %
        %   Example:
        %       Problem = PROBLEM.Current()
        
            persistent Problem;
            if nargin > 0
                Problem = obj;
            end
            if nargout > 0
                obj = Problem;
            end
        end
    end
	methods(Access = protected, Sealed)
        function varargout = ParameterSet(obj,varargin)
        %ParameterSet - Obtain the parameters of the problem.
        %
        %   [p1,p2,...] = obj.ParameterSet(v1,v2,...) sets the values of
        %   parameters p1, p2, ..., where each parameter is set to the
        %   value given in obj.parameter if obj.parameter is specified, and
        %   set to the value given in v1, v2, ... otherwise.
        %
        %   Example:
        %       [p1,p2,p3] = obj.ParameterSet(1,2,3)

            varargout = varargin;
            specified = ~cellfun(@isempty,obj.parameter);
            varargout(specified) = obj.parameter(specified);
        end
    end
end