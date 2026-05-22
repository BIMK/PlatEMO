classdef SOLUTION_SparseEMT < handle
%SOLUTION_SparseEMT - The class of a solution of SparseEMT.
%
%   This is the class of a solution. An object of SOLUTION stores all the
%   properties including decision variables, objective values, constraint
%   violations, and additional properties of a solution.
%
% SOLUTION_SparseEMT properties:
%   dec         <vector>	real decision variables of the solution 
%   mask        <vector>	binary decision variables of the solution
%   dec_high    <vector>  	high real ecision variables of the solution 
%   mask_high   <vector>  	high binary decision variables of the solution  
%   obj         <vector>   	objective values of the solution
%   con         <vector>  	constraint violations of the solution      
%   rank        <vector>  	rank of the solution
%   skill_factor<vector>  	skill factor of the individual
%   add         <vector>   	additional properties of the solution
%
% SOLUTION methods:
%   SOLUTION	<private>  	the constructor setting all the properties of 
%                         	solutions
%   decs        <public>  	get the matrix of real decision variables of
%                          	multiple solutions
%   masks       <public>  	get the matrix of binary decision variables of
%                           multiple solutions
%   objs        <public>   	get the matrix of objective values of multiple 
%                          	solutions
%   cons        <public>  	get the matrix of constraint violations of
%                          	multiple solutions
%   adds        <public>  	get the matrix of additional properties of
%                         	multiple solutions
%   best        <public>   	get the feasible and nondominated solutions
%                          	among multiple solutions

    properties(SetAccess = public)
        dec;        % Real decision variables of the solution
        mask;       % Binary decision variables of the solution
        dec_high;   % High real decision variables of the solution
        mask_high;  % High binary decision variables of the solution       
        obj;        % Objective values of the solution
        con;        % Constraint violations of the solution    
        rank;       % rank of the solution
        skill_factor;        % skill factor of the solution
    end
    methods
        %% Constructor
        function obj = SOLUTION_SparseEMT(Decs,Masks,Decs_high,Masks_high,Tasks,skill_factor,Problem)
        %SOLUTION_SparseEMT - Constructor of SOLUTION_SparseEMT class
        
            if nargin > 0
                % Create new objects
                obj(1,size(Decs,1)) = SOLUTION_SparseEMT;
                Upper = repmat(Tasks(skill_factor).upper,length(obj),1);
                Lower = repmat(Tasks(skill_factor).lower,length(obj),1);
                Decs  = max(min(Decs,Upper),Lower);

                Upper_high = repmat(Tasks(end).upper,length(obj),1);
                Lower_high = repmat(Tasks(end).lower,length(obj),1);
                Decs_high  = max(min(Decs_high,Upper_high),Lower_high);
                              
                % Assign the decision variables, objective values,
                % constraint violations, and additional properties
                for i = 1 : length(obj)
                    obj(i).con  = 0;
                    obj(i).dec  = Decs(i,:);
                    obj(i).mask = Masks(i,:);
                    obj(i).dec_high  = Decs_high(i,:);
                    obj(i).mask_high = Masks_high(i,:);
                    obj(i).obj  = Problem.Evaluation(Decs_high(i,:).*Masks_high(i,:)).obj;
                    obj(i).skill_factor = skill_factor;
                end
            end
        end
        %% Get the matrix of decision variables of the population
        function value = decs(obj)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of 
        %   the population obj, where obj is an array of INDIVIDUAL objects.
            Decs = [];
            for i = 1:length(obj)
                Decs(i,:) = obj(i).dec_high.*obj(i).mask_high;
            end
            value = cat(1,Decs);
        end
        %% Get the matrix of real decision variables of the population
        function value = decss(obj)
        %decss - Get the matrix of real decision variables of the population.
        %
        %   A = obj.decss returns the matrix of real decision variables of 
        %   the population obj, where obj is an array of INDIVIDUAL objects.
            Decs = [];
            for i = 1:length(obj)
                Decs(i,:) = obj(i).dec;
            end
            value = cat(1,Decs);
        end
        %% Get the matrix of binary decision variables of the population
        function value = masks(obj)
        %masks - Get the matrix of binary decision variables of the population.
        %
        %   A = obj.masks returns the matrix of binary decision variables of 
        %   the population obj, where obj is an array of INDIVIDUAL objects.
            Masks = [];
            for i = 1:length(obj)
                Masks(i,:) = obj(i).mask;
            end
            value = cat(1,Masks);
        end
        %% Get the matrix of high-dimensional real decision variables of the population
        function value = decss_high(obj)
        %decss - Get the matrix of real decision variables of the population.
        %
        %   A = obj.decss returns the matrix of real decision variables of 
        %   the population obj, where obj is an array of INDIVIDUAL objects.
            Decs = [];
            for i = 1:length(obj)
                Decs(i,:) = obj(i).dec_high;
            end
            value = cat(1,Decs);
        end
        %% Get the matrix of high-dimensional binary decision variables of the population
        function value = masks_high(obj)
        %masks - Get the matrix of binary decision variables of the population.
        %
        %   A = obj.masks returns the matrix of binary decision variables of 
        %   the population obj, where obj is an array of INDIVIDUAL objects.
            Masks = [];
            for i = 1:length(obj)
                Masks(i,:) = obj(i).mask_high;
            end
            value = cat(1,Masks);
        end
        %% Get the matrix of objective values of the population
        function value = objs(obj)
        %objs - Get the matrix of objective values of the population.
        %
        %   A = obj.objs returns the matrix of objective values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.obj);
        end
        %% Get the matrix of constraint violations of the population
        function value = cons(obj)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.con);
        end
        %% Get the matrix of additional properties of the population
        function value = skill_factors(obj,AddProper)
        %skill_factors - Get the matrix of skill factors of the population.
        %   A = obj.skill_factors returns the matrix of skill factors of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            for i = 1 : length(obj)
                if isempty(obj(i).skill_factor)
                    obj(i).skill_factor = AddProper(i,:);
                end
            end
            value = cat(1,obj.skill_factor);
        end
        function value = ranks(obj,AddProper)
        %ranks - Get the matrix of ranks of the population.
        %   A = obj.ranks returns the matrix of ranks of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            for i = 1 : length(obj)
                if isempty(obj(i).rank)
                    obj(i).rank = AddProper(i,:);
                end
            end
            value = cat(1,obj.rank);
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