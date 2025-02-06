classdef RBFN_P
    
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Haoran Gu

    properties
        Problem;
        index;
        x;
        trainsample;
        vec;
        d;
        lam;
    end
    methods
        %% Constructor
        function obj = RBFN_P(Problem)
            if nargin == 1
                obj(1, Problem.N) = RBFN_P;
                for i = 1 : length(obj)
                    obj(i).index   = i;
                    obj(i).Problem = Problem;
                    obj(i).x       = [];
                    obj(i).trainsample;
                    obj(i).vec;
                    obj(i).d;
                    obj(i).lam;
                end
            end
        end
        %% Model-construction
        function obj = ModelConstruction(obj, A, B_i, W, Z,center)
            indices = 1 : length(A);
            for i = 1 : length(A)
                obj.x(i, :) = A(i).dec;
            end
            g_data1 = max(abs(A(indices).objs - repmat(Z, length(indices), 1)) .* W(center, :), [], 2);
            uniformed_xdata = zeros(length(A), obj.Problem.D);
            for i = 1 : length(A)
                uniformed_xdata(i, :) = obj.UniformInput(obj.x(i, :));
            end
            covx            = cov(uniformed_xdata);
            [vec_pop,~,~]   = pcacov(covx);
            uniformed_xdata = uniformed_xdata* vec_pop;
            obj.vec         = vec_pop;
            
            obj.d = obj.Problem.D/2;
            obj.trainsample = uniformed_xdata(:,1:obj.d);
            pair            = pdist2(obj.trainsample, obj.trainsample);
            D_max           = max(max(pair, [], 2));
            spread          = D_max * (obj.Problem.D * obj.Problem.N) ^ (-1 / obj.Problem.D);
            RBFModel        = newrbe(transpose(obj.trainsample), transpose(g_data1), spread);
            obj.lam         = RBFModel;
        end
        function predicted_value = PredictClass(obj, x)
            uniformed_xdata = zeros(size(obj.x,1), obj.Problem.D);
            uniformed_x     = zeros(size(x,1), obj.Problem.D);
            for i = 1 : size(obj.x,1)
                uniformed_xdata(i, :) = obj.UniformInput(obj.x(i, :));
            end
            for j = 1 : size(x,1)
                uniformed_x(j,:) = obj.UniformInput(x(j,:));
            end
            uniformed_x     = uniformed_x * obj.vec;
            test_x          = uniformed_x(:,1:obj.d);
            predicted_value = transpose(sim(obj.lam, transpose(test_x)));
        end
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