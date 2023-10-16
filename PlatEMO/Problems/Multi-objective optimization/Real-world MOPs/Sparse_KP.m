classdef Sparse_KP < PROBLEM
% <multi/many> <binary> <large/none>
% The sparse multi-objective knapsack problem

%------------------------------- Reference --------------------------------
% Y. Su, Z. Jin, Y. Tian, X. Zhang, and K. C. Tan, Comparing the
% performance of evolutionary algorithms for sparse multi-objective
% optimization via a comprehensive indicator, IEEE Computational
% Intelligence Magazine, 2022, 17(3): 34-53.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        P;	% Profit of each item according to each knapsack
        W;  % Weight of each item according to each knapsack
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 250; end
            obj.encoding = 4 + zeros(1,obj.D);
            % Randomly generate profits and weights
            file = sprintf('Dataset_KP-M%d-D%d.mat',obj.M,obj.D);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'P','W');
            else
                P = randi([10,100],obj.M,obj.D);
                W = randi([10,100],obj.M,obj.D);
                save(file,'P','W');
            end
            obj.P = P;
            obj.W = W;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = repmat(sum(obj.P,2)',size(PopDec,1),1) - PopDec*obj.P';
            PopObj = PopObj + 10*repmat(max(max(0,PopDec*obj.W'-repmat(0.1*sum(obj.W,2)',size(PopDec,1),1)),[],2),1,size(PopObj,2));
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = sum(obj.P,2)';
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            Draw(Population.decs*obj.P',{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
        end
    end
end
