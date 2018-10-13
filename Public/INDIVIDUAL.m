classdef INDIVIDUAL < handle
%INDIVIDUAL - The class of an individual
%
%   This is the class of an individual. An object of INDIVIDUAL stores all
%   the properties including decision variables, objective values,
%   constraint values and additional properties of an individual.
%
% INDIVIDUAL properties:
%   dec         <read-only>     decision variables of the individual
%   obj         <read-only>     objective values of the individual
%   con         <read-only>     constraint values of the individual
%   add         <read-only>     additional properties of the individual
%
% INDIVIDUAL methods:
%   INDIVIDUAL	<public>        the constructor, all the properties will be
%                               set when the object is creating
%   decs        <public>      	get the matrix of decision variables of the
%                               population
%   objs        <public>        get the matrix of objective values of the
%                               population
%   cons        <public>        get the matrix of constraint values of the
%                               population
%   adds        <public>        get the matrix of additional property of
%                               the population

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        dec;        % Decision variables of the individual
        obj;        % Objective values of the individual
        con;        % Constraint values of the individual
        add;        % Additional properties of the individual
    end
    methods
        %% Constructor
        function obj = INDIVIDUAL(Decs,addValue)
        %INDIVIDUAL - Constructor of INDIVIDUAL class
        %
        %   H = INDIVIDUAL(Dec) creates an array of individuals (i.e. a
        %   population), where Dec is the matrix of decision variables of
        %   the population. The objective values and the constraint values
        %   are automatically calculated by the test problem functions.
        %   After creating the objects, the number of evaluations will be
        %   increased by length(H).
        %
        %   H = INDIVIDUAL(Dec,AddProper) creates the population with
        %   additional properties stored in AddProper. The name of the
        %   additional property is the same to the function name of the
        %   caller.
        %
        %   Example:
        %       H = INDIVIDUAL(rand(100,3))
        %       H = INDIVIDUAL(rand(100,10),randn(100,3))
        
            if nargin > 0
                % Create new objects
                obj(1,size(Decs,1)) = INDIVIDUAL;
                Global = INDIVIDUAL.register;
                if ~isempty(Global)
                    % Set the infeasible decision variables to boundary values
                    if ~isempty(Global.lower) && ~isempty(Global.upper)
                        Lower = repmat(Global.lower,length(obj),1);
                        Upper = repmat(Global.upper,length(obj),1);
                        Decs  = max(min(Decs,Upper),Lower);
                    end
                    % Calculte the objective and constraint values
                    [Decs,Objs,Cons] = Global.problem('value',Global,Decs);
                    if isempty(Cons)
                        Cons = zeros(length(obj),1);
                    end
                    % Assign the decision variables, objective values,
                    % constraint values and additional property values
                    if nargin < 2
                        for i = 1 : length(obj)
                            obj(i).dec = Decs(i,:);
                            obj(i).obj = Objs(i,:);
                            obj(i).con = Cons(i,:);
                        end
                    else
                        CallStack = dbstack();
                        Field     = CallStack(2).name;
                        for i = 1 : length(obj)
                            obj(i).dec         = Decs(i,:);
                            obj(i).obj         = Objs(i,:);
                            obj(i).con         = Cons(i,:);
                            obj(i).add.(Field) = addValue(i,:);
                        end
                    end
                    % Update the number of evaluated individuals
                    Global.evaluated = Global.evaluated + length(obj);
                else
                    % Assign the decision variables
                    for i = 1 : length(obj)
                        obj(i).dec = Decs(i,:);
                    end
                end
            end
        end
        %% Get the matrix of decision variables of the population
        function value = decs(obj)
        %decs - Get the matrix of decision variables of the population
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.dec);
        end
        %% Get the matrix of objective values of the population
        function value = objs(obj)
        %objs - Get the matrix of objective values of the population
        %
        %   A = obj.objs returns the matrix of objective values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.obj);
        end
        %% Get the matrix of constraint values of the population
        function value = cons(obj)
        %cons - Get the matrix of constraint values of the population
        %
        %   A = obj.cons returns the matrix of constraint values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.con);
        end
        %% Get the matrix of additional property of the population
        function value = adds(obj,addValue)
        %adds - Get the matrix of additional property values of the population
        %
        %   A = obj.adds(AddProper) returns the matrix of the values of the
        %   additional property of the INDIVIDUAL objects obj. The name of
        %   the additional property is same to the function name of the
        %   caller, that is, the values of one additional property of the
        %   individuals can only be obtained by the function which created
        %   them. If any individual in obj does not contain the specified
        %   additional property, assign it a default value specified in
        %   AddProper.
        
            CallStack = dbstack();
            Field     = CallStack(2).name;
            value     = zeros(length(obj),size(addValue,2));
            for i = 1 : length(obj)
                if ~isfield(obj(i).add,Field)
                    obj(i).add.(Field) = addValue(i,:);
                end
                value(i,:) = obj(i).add.(Field);
            end
        end
    end
    methods(Static, Access = ?GLOBAL)
        %% Register a GLOBAL object with INDIVIDUAL class
        function obj = register(obj)
            persistent Global;
            if nargin > 0
                Global = obj;
            else
                obj = Global;
            end
        end
    end
end