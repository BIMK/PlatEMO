classdef SVM
% Definition of SVM model and fuctions used in MCEA/D

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Masaya Nakata

    properties 
        Problem           % Problem instance defined by PlatEMO
        index             % Index of sub-problem
        x                 % SVM input
        label             % Class of training input
        mdl               % SVM model
        C                 % SVM parameter C
        gamma             % SVM parameter gamma
    end

    methods
        %% Constructor
        function obj = SVM(Problem)
            if nargin == 1
                obj(1, Problem.N) = SVM;
                for i = 1 : length(obj)
                    obj(i).index        = i;
                    obj(i).Problem      = Problem;
                    obj(i).C            = 1.0;
                    obj(i).gamma        = 1.0;
                    obj(i).x            = [];
                    obj(i).label        = [];
                    obj(i).mdl;
                end
            end
        end

        %% Model-construction
        function obj = ModelConstruction(obj, A, B_i, W, Z)
            % Initialization
            indices = [1 : length(A)];
            for i = 1 : length(A)
                obj.x(i, :)     = A(i).dec;
                obj.label(i, 1) = -1;
            end
            
            % Get the set of current best solutions of neighbor sub-problems
            C_i = [];
            for i = 1 : length(B_i)
                % Calculate scalarization funcion values
                g_data            = max(abs(A(indices).objs - repmat(Z, length(indices), 1)) .* W(B_i(i), :), [], 2);
                
                % Get the set of current best solution avoiding duplicative selection
                [~, sorted_index] = sort(g_data);
                for j = 1 : length(sorted_index)
                    if ~ismember(sorted_index(j), C_i)
                        C_i = [C_i, sorted_index(j)];
                        obj.label(sorted_index(j), 1) = 1;
                        break
                    end
                end
            end

            % Train SVM 
            uniformed_xdata = zeros(length(obj.label), obj.Problem.D);
            for i = 1 : length(obj.label)
                uniformed_xdata(i, :) = obj.UniformInput(obj.x(i, :));
            end
            sigma           = sqrt(1 / (2 * obj.gamma));
            obj.mdl         = fitcsvm(uniformed_xdata, obj.label, 'BoxConstraint', obj.C, 'KernelScale', sigma, 'KernelFunction', 'rbf');
        end

        %% Predict the class of input and get the decision score function value
        function [predicted_class, score] = PredictClass(obj, x)
            % Uniform the input
            uniformed_x = obj.UniformInput(x);
            
            % Predict the class of input
            [predicted_class, score_list] = obj.mdl.predict(uniformed_x);
            
            % Return the decision score function value
            score = score_list(2);
        end

        %% Uniform the input
        function uniformed_x = UniformInput(obj, x)
            uniformed_x = ones(1, obj.Problem.D);
            for i = 1 : obj.Problem.D
                x_min = obj.Problem.lower(i);
                x_max = obj.Problem.upper(i);
                uniformed_x(i) = (x(i) - x_min) / (x_max - x_min);
            end
        end
        
    end
end