classdef SOLUTION < handle
%SOLUTION - The class of a solution.
%
%   This is the class of a solution. An object of SOLUTION stores all the
%   properties including decision variables, objective values, constraint
%   violations, and additional properties of a solution.
%
% SOLUTION properties:
%   dec         <vector>	decision variables of the solution
%   obj         <vector>   	objective values of the solution
%   con         <vector>  	constraint violations of the solution
%   add         <vector>   	additional properties of the solution
%
% SOLUTION methods:
%   SOLUTION	<private>  	the constructor setting all the properties of 
%                         	solutions
%   decs        <public>  	get the matrix of decision variables of
%                          	multiple solutions
%   objs        <public>   	get the matrix of objective values of multiple 
%                          	solutions
%   cons        <public>  	get the matrix of constraint violations of
%                          	multiple solutions
%   adds        <public>  	get the matrix of additional properties of
%                         	multiple solutions
%   best        <public>   	get the feasible and nondominated solutions
%                          	among multiple solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
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
    end
    properties
        add;        % Additional properties of the solution
    end
    methods
        function obj = SOLUTION(PopDec,PopObj,PopCon,PopAdd)
        %SOLUTION - The constructor of SOLUTION.
        %
        %   P = SOLUTION(Dec,Obj,Con) creates an array of SOLUTION objects
        %   with decision variables of Dec, objective values of Obj, and
        %   constraint violations of Con.
        %
        %   P = SOLUTION(Dec,Obj,Con,Add) also sets the additional
        %   properties (e.g., velocity) of solutions.
        %
        %   Dec, Obj, Con, Add are matrices, where each row denotes a
        %   solution and each column denotes a dimension of the variables,
        %   objectives, constraints, or additional properties.
        %
        %   Example:
        %       Population = SOLUTION(PopDec,PopObj,PopCon)
        %       Population = SOLUTION(PopDec,PopObj,PopCon,PopVel)
        
            if nargin > 0
                obj(1,size(PopDec,1)) = SOLUTION;
                for i = 1 : length(obj)
                    obj(i).dec = PopDec(i,:);
                    obj(i).obj = PopObj(i,:);
                    obj(i).con = PopCon(i,:);
                end
                if nargin > 3
                    for i = 1 : length(obj)
                        obj(i).add = PopAdd(i,:);
                    end
                end
            end
        end
    end
    methods
        function value = decs(obj)
        %decs - Get the matrix of decision variables of multiple solutions.
        %
        %   Dec = obj.decs returns the matrix of decision variables of
        %   multiple solutions obj.
        
            value = cat(1,obj.dec);
        end
        function value = objs(obj)
        %objs - Get the matrix of objective values of multiple solutions.
        %
        %   Obj = obj.objs returns the matrix of objective values of
        %   multiple solutions obj.
        
            value = cat(1,obj.obj);
        end
        function value = cons(obj)
        %cons - Get the matrix of constraint violations of multiple solutions.
        %
        %   Con = obj.cons returns the matrix of constraint violations of
        %   multiple solutions obj.
        
            value = cat(1,obj.con);
        end
        function value = adds(obj,Add)
        %adds - Get or set the matrix of additional properties of multiple solutions.
        %
        %   Add = obj.adds(Add) returns the matrix of additional properties
        %   of multiple solutions obj. If any solution in obj does not
        %   contain an additional property, it will be set to the default
        %   value specified in Add.

            for i = 1 : length(obj)
                if isempty(obj(i).add)
                    obj(i).add = Add(i,:);
                end
            end
            value = cat(1,obj.add);
        end
        function P = best(obj)
        %best - Get the best solutions among multiple solutions.
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