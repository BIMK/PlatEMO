classdef PROBLEM < handle & matlab.mixin.Heterogeneous
%PROBLEM - The superclass of problems.
%
%   This is the superclass of problems. An object of PROBLEM stores all the
%   settings of the problem.
%
% PROBLEM properties:
%   N               <scalar>    population size
%   M               <scalar>    number of objectives
%   D               <scalar>    number of decision variables
%   maxFE           <scalar>    maximum number of function evaluations
%   FE              <scalar>    number of consumed function evaluations
%   maxRuntime      <scalar>    maximum runtime (in second)
%   encoding        <vector>    encoding scheme of each decision variable (1.real 2.integer 3.label 4.binary 5.permutation)
%   lower           <vector>    lower bound of each decision variable
%   upper           <vector>   	upper bound of each decision variable
%   optimum         <matrix>    optimal objective values of the problem
%   PF              <matrix>    image of the Pareto front
%   parameter       <any>       other parameters of the problem
%
% PROBLEM methods:
%   PROBLEM         <protected> the constructor setting all the properties specified by user
%   Setting         <public>    default settings of the problem
%   Initialization 	<public>    generate initial solutions
%   Evaluation      <public>    evaluate solutions
%   CalDec          <public>    repair invalid solutions
%   CalObj          <public>    calculate the objective values of solutions
%   CalCon          <public>    calculate the constraint violations of solutions
%   CalObjGrad      <public>    calculate the gradients of objectives
%   CalConGrad      <public>    calculate the gradients of constraints
%   GetOptimum      <public>    generate the optimal objective values of the problem
%   GetPF          	<public>    generate the image of the Pareto front
%   CalMetric       <public>    calculate the metric value of a population
%   DrawDec         <public>    display solutions in the decision space
%   DrawObj         <public>    display solutions in the objective space
%   ParameterSet	<protected>	obtain the parameter settings of the problem

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        N          = 100;      	% Population size
        maxFE      = 10000;     % Maximum number of function evaluations
        FE         = 0;        	% Number of consumed function evaluations
    end
    properties(SetAccess = protected)
        M;                    	% Number of objectives
        D;                     	% Number of decision variables
        maxRuntime = inf;      	% maximum runtime (in second)
        encoding   = 1;        	% Encoding scheme of each decision variable (1.real 2.integer 3.label 4.binary 5.permutation)
        lower      = 0;     	% Lower bound of each decision variable
        upper      = 1;        	% Upper bound of each decision variable
        optimum;              	% Optimal values of the problem
        PF;                   	% Image of Pareto front
        parameter  = {};       	% Other parameters of the problem
    end
    methods(Access = protected)
        function obj = PROBLEM(varargin)
        %PROBLEM - The constructor of PROBLEM.
        %
        %   Problem = proName('Name',Value,'Name',Value,...) generates an
        %   object with the properties specified by the inputs. proName is
        %   a subclass of PROBLEM, while PROBLEM cannot be instantiated
        %   directly.
        %
        %   If proName is UserProblem, all the properties can be specified
        %   to define the details of the problem. Otherwise, only the
        %   properties N, M, D, maxFE, maxRuntime can be defined. The
        %   properties M, D, encoding, lower, upper may be automatically
        %   revised after specification.
        %
        %   Example:
        %       Problem = UserProblem('objFcn',@(x)sum(x,2))
        %       Problem = DTLZ2('M',5,'D',10)

            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'N','M','D','maxFE','maxRuntime','parameter'}))
                obj.(varargin{i}) = varargin{i+1};
            end
            obj.Setting();
            obj.optimum  = obj.GetOptimum(10000);
            obj.PF       = obj.GetPF();
        end
    end
    methods
        function Setting(obj)
        %Setting - Default settings of the problem.
        %
        %   This function is expected to be implemented in each subclass of
        %   PROBLEM, which is usually called by the constructor.
        end
        function Population = Initialization(obj,N)
        %Initialization - Generate multiple initial solutions.
        %
        %   P = obj.Initialization() randomly generates the decision
        %   variables of obj.N solutions and returns the SOLUTION objects.
        %
        %   P = obj.Initialization(N) generates N solutions.
        %
        %   This function is usually called at the beginning of algorithms.
        %
        %   Example:
        %       Population = Problem.Initialization()
        
            if nargin < 2
            	N = obj.N;
            end
            PopDec = zeros(N,obj.D);
            Type   = arrayfun(@(i)find(obj.encoding==i),1:5,'UniformOutput',false);
            if ~isempty(Type{1})        % Real variables
                PopDec(:,Type{1}) = unifrnd(repmat(obj.lower(Type{1}),N,1),repmat(obj.upper(Type{1}),N,1));
            end
            if ~isempty(Type{2})        % Integer variables
                PopDec(:,Type{2}) = round(unifrnd(repmat(obj.lower(Type{2}),N,1),repmat(obj.upper(Type{2}),N,1)));
            end
            if ~isempty(Type{3})        % Label variables
                PopDec(:,Type{3}) = round(unifrnd(repmat(obj.lower(Type{3}),N,1),repmat(obj.upper(Type{3}),N,1)));
            end
            if ~isempty(Type{4})        % Binary variables
                PopDec(:,Type{4}) = logical(randi([0,1],N,length(Type{4})));
            end
            if ~isempty(Type{5})        % Permutation variables
                [~,PopDec(:,Type{5})] = sort(rand(N,length(Type{5})),2);
            end
            Population = obj.Evaluation(PopDec);
        end
        function Population = Evaluation(obj,varargin)
        %Evaluation - Evaluate multiple solutions.
        %
        %   P = obj.Evaluation(Dec) returns the SOLUTION objects based on
        %   the decision variables Dec. The objective values and constraint
        %   violations of the solutions are calculated automatically, and
        %   obj.FE is increased accordingly.
        %
        %   P = obj.Evaluation(Dec,Add) also sets the additional properties
        %   (e.g., velocity) of solutions.
        %
        %   This function is usually called after generating new solutions.
        %
        %   Example:
        %       Population = Problem.Evaluation(PopDec)
        %       Population = Problem.Evaluation(PopDec,PopVel)
        
            PopDec     = obj.CalDec(varargin{1});
            PopObj     = obj.CalObj(PopDec);
            PopCon     = obj.CalCon(PopDec);
            Population = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        function PopDec = CalDec(obj,PopDec)
        %CalDec - Repair multiple invalid solutions.
        %
        %   Dec = obj.CalDec(Dec) repairs the invalid (not infeasible)
        %   decision variables in Dec.
        %
        %   An invalid solution indicates that it is out of the decision
        %   space, while an infeasible solution indicates that it does not
        %   satisfy all the constraints.
        %
        %   This function is usually called by PROBLEM.Evaluation.
        %
        %   Example:
        %       PopDec = Problem.CalDec(PopDec)

            Type  = arrayfun(@(i)find(obj.encoding==i),1:5,'UniformOutput',false);
            index = [Type{1:3}];
            if ~isempty(index)
                PopDec(:,index) = max(min(PopDec(:,index),repmat(obj.upper(index),size(PopDec,1),1)),repmat(obj.lower(index),size(PopDec,1),1));
            end
            index = [Type{2:5}];
            if ~isempty(index)
                PopDec(:,index) = round(PopDec(:,index));
            end
        end
        function PopObj = CalObj(obj,PopDec)
        %CalObj - Calculate the objective values of multiple solutions.
        %
        %   Obj = obj.CalObj(Dec) returns the objective values of Dec.
        %
        %   This function is usually called by PROBLEM.Evaluation.
        %
        %   Example:
        %       PopObj = Problem.CalObj(PopDec)

            PopObj = zeros(size(PopDec,1),1);
        end
        function PopCon = CalCon(obj,PopDec)
        %CalCon - Calculate the constraint violations of multiple solutions.
        %
        %   Con = obj.CalCon(Dec) returns the constraint violations of Dec.
        %
        %   This function is usually called by PROBLEM.Evaluation.
        %
        %   Example:
        %       PopCon = Problem.CalCon(PopDec)
        
            PopCon = zeros(size(PopDec,1),1);
        end
        function ObjGrad = CalObjGrad(obj,Dec)
        %CalObjGrad - Calculate the gradients of objectives of a solution.
        %
        %   Grad = obj.CalObjGrad(Dec) returns the gradients of objectives
        %   of Dec, i.e., a Jacobian matrix.
        %
        %   This function is usually called by gradient-based algorithms.
        %
        %   Example:
        %       ObjGrad = Problem.CalObjGrad(Dec)

            Dec(Dec==0) = 1e-12;
            X           = repmat(Dec,length(Dec),1).*(1+eye(length(Dec))*1e-6);
            ObjGrad     = (obj.CalObj(X)-repmat(obj.CalObj(Dec),size(X,1),1))'./Dec./1e-6;
        end
        function ConGrad = CalConGrad(obj,Dec)
        %CalConGrad - Calculate the gradients of constraints of a solution.
        %
        %   Grad = obj.CalConGrad(Dec) returns the gradients of constraints
        %   of Dec, i.e., a Jacobian matrix.
        %
        %   This function is usually called by gradient-based algorithms.
        %
        %   Example:
        %       ConGrad = Problem.CalConGrad(Dec)
        
            Dec(Dec==0) = 1e-12;
            X           = repmat(Dec,length(Dec),1).*(1+eye(length(Dec))*1e-6);
            ConGrad     = (obj.CalCon(X)-repmat(obj.CalCon(Dec),size(X,1),1))'./Dec./1e-6;
        end
        function R = GetOptimum(obj,N)
        %GetOptimum - Generate the optimums of the problem.
        %
        %   R = obj.GetOptimum(N) returns N optimums of the problem for
        %   metric calculation.
        %
        %   For single-objective optimization problems, an optimum can be
        %   the minimum objective value of the problem.
        %
        %   For multi-objective optimization problems, an optimum can be a
        %   point on the Pareto front; if the Pareto front is unknown, an
        %   optimum can be a reference point for hypervolume calculation.
        %
        %   This function is usually called by the constructor.
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
        %   For single-objective optimization problems, this function is
        %   useless.
        %
        %   For bi-objective optimization problems, the image should be a
        %   one-dimensional curve.
        %
        %   For tri-objective optimization problems, the image should be a
        %   two-dimensional surface.
        %
        %   For constrained bi-objective optimization problems, the image
        %   can be the feasible regions.
        %
        %   This function is usually called by the constructor.
        %
        %   Example:
        %       R = Problem.GetPF()
        
            R = [];
        end
        function score = CalMetric(obj,metName,Population)
        %CalMetric - calculate the metric value of a population
        %
        %   value = obj.CalMetric(Met,P) returns the metric value of a
        %   population P, where Met is a string denoting the name of a
        %   metric function.
        %
        %   Example:
        %       value = Problem.CalMetric('HV',Population);
        
            score = feval(metName,Population,obj.optimum);
        end
        function DrawDec(obj,Population)
        %DrawDec - Display a population in the decision space.
        %
        %   obj.DrawDec(P) displays the decision variables of population P.
        %
        %   This function is usually called by the GUI.
        %
        %   Example:
        %       Problem.DrawDec(Population)
        
            if all(obj.encoding==4)
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
        %	This function is usually called by the GUI.
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
	methods(Access = protected, Sealed)
        function varargout = ParameterSet(obj,varargin)
        %ParameterSet - Obtain the parameters of the problem.
        %
        %   [p1,p2,...] = obj.ParameterSet(v1,v2,...) sets the values of
        %   parameters p1, p2, ..., where each parameter is set to the
        %   value given in obj.parameter if obj.parameter is specified, and
        %   set to the value given in v1, v2, ... otherwise.
        %
        %   This function is usually called by PROBLEM.Setting.
        %
        %   Example:
        %       [p1,p2,p3] = obj.ParameterSet(1,2,3)

            varargout = varargin;
            specified = ~cellfun(@isempty,obj.parameter);
            varargout(specified) = obj.parameter(specified);
        end
    end
end