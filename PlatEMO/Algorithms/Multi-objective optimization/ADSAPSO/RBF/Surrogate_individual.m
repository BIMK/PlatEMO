classdef Surrogate_individual
%Surrogate_individual - The class of an surrogate individual.

% This function is written by Jianqing Lin

    properties
        dec;        % Decision variables of the individual
        obj;        % Objective values of the individual
        add;        % Additional properties of the individual
    end
    methods
        %% Constructor
        function obj = Surrogate_individual(Decs,Objs,AddProper)
            if nargin > 0
                % Create new objects
                obj(1,size(Decs,1)) = Surrogate_individual;
                % Assign the decision variables, objective values and additional properties
                for i = 1 : length(obj)
                    obj(i).dec = Decs(i,:);
                    obj(i).obj = Objs(i,:);
                end
                if nargin > 2
                    for i = 1 : length(obj)
                        obj(i).add = AddProper(i,:);
                    end
                end
            end
        end
        %% Get the matrix of decision variables of the population
        function value = decs(obj)
        %decs - Get the matrix of decision variables of the population.
        
            value = cat(1,obj.dec);
        end
        %% Get the matrix of objective values of the population
        function value = objs(obj)
        %objs - Get the matrix of objective values of the population.
        
            value = cat(1,obj.obj);
        end
        %% Get the matrix of additional properties of the population
        function value = adds(obj,AddProper)
        %adds - Get the matrix of additional properties of the population.

            for i = 1 : length(obj)
                if isempty(obj(i).add)
                    obj(i).add = AddProper(i,:);
                end
            end
            value = cat(1,obj.add);
        end
    end
end