classdef SOLUTION < handle
%SOLUTION - The class of a solution.
%
%   This is the class of a solution. An object of SOLUTION stores all the
%   properties including decision variables, objective values, constraint
%   violations, and additional properties of a solution.
%
% SOLUTION properties:
%   dec         <read-only>     decision variables of the solution
%   obj         <read-only>     objective values of the solution
%   con         <read-only>     constraint violations of the solution
%   add         <read-only>     additional properties of the solution
%
% SOLUTION methods:
%   SOLUTION	<public>        the constructor, which sets all the
%                               properties of the solution
%   decs        <public>      	get the matrix of decision variables of
%                               multiple solutions
%   objs        <public>        get the matrix of objective values of
%                               multiple solutions
%   cons        <public>        get the matrix of constraint violations of
%                               multiple solutions
%   adds        <public>        get the matrix of additional properties of
%                               multiple solutions
%   best        <public>        get the feasible and nondominated solutions
%                               among multiple solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        dec;        % Decision variables of the solution
        obj;        % Objective values of the solution
        con;        % Constraint violations of the solution
        add;        % Additional properties of the solution
    end
    methods
        function obj = SOLUTION(PopDec,AddPro)
        %SOLUTION - The constructor of SOLUTION.
        %
        %   P = SOLUTION(Dec) creates an array of SOLUTION objects with
        %   decision variables of Dec, where the objective values and
        %   constraint violations are calculated automatically.
        %
        %   P = SOLUTION(Dec,AddPro) also sets the additional properties
        %   (e.g., velocity) of solutions to the values of AddPro.
        %
        %   Dec and AddPro are matrices, where each row denotes a solution
        %   and each column denotes a dimension of the decision variables
        %   or additional properties.
        %
        %   Example:
        %       Population = SOLUTION(PopDec)
        
            if nargin > 0
                obj(1,size(PopDec,1)) = SOLUTION;
                Problem = PROBLEM.Current();
                PopDec  = Problem.CalDec(PopDec);
                PopObj  = Problem.CalObj(PopDec);
                PopCon  = Problem.CalCon(PopDec);
                for i = 1 : length(obj)
                    obj(i).dec = PopDec(i,:);
                    obj(i).obj = PopObj(i,:);
                    obj(i).con = PopCon(i,:);
                end
                if nargin > 1
                    for i = 1 : length(obj)
                        obj(i).add = AddPro(i,:);
                    end
                end
                Problem.FE = Problem.FE + length(obj);
            end
        end
        function value = decs(obj)
        %decs - Get the matrix of decision variables of a population.
        %
        %   Dec = obj.decs returns the matrix of decision variables of
        %   multiple solutions obj.
        
            value = cat(1,obj.dec);
        end
        function value = objs(obj)
        %objs - Get the matrix of objective values of a population.
        %
        %   Obj = obj.objs returns the matrix of objective values of
        %   multiple solutions obj.
        
            value = cat(1,obj.obj);
        end
        function value = cons(obj)
        %cons - Get the matrix of constraint violations of a population.
        %
        %   Con = obj.cons returns the matrix of constraint violations of
        %   multiple solutions obj.
        
            value = cat(1,obj.con);
        end
        function value = adds(obj,AddPro)
        %adds - Get the matrix of additional properties of a population.
        %
        %   Add = obj.adds(AddPro) returns the matrix of additional
        %   properties of multiple solutions obj. If any solution in obj
        %   does not contain an additional property, it will be set to the
        %   default value specified in AddPro.

            for i = 1 : length(obj)
                if isempty(obj(i).add)
                    obj(i).add = AddPro(i,:);
                end
            end
            value = cat(1,obj.add);
        end
        function P = best(obj)
        %best - Get the best solutions in a population.
        %
        %   P = obj.best returns the feasible and non-dominated solutions
        %   among multiple solutions obj. If the solutions have a single
        %   objective, the feasible solution with minimum objective value
        %   is returned.
        
            Feasible = find(all(obj.cons<=0,2));
            if isempty(Feasible)
                Best = [];
            elseif length(obj(1).obj) > 1
                Best = NDSort(obj(Feasible).objs,1) == 1;
            else
                [~,Best] = min(obj(Feasible).objs);
            end
            P = obj(Feasible(Best));
        end
    end
end