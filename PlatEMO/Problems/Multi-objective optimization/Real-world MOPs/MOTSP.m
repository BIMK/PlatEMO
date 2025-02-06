classdef MOTSP < PROBLEM
% <2007> <multi/many> <permutation> <large/none>
% The multi-objective traveling salesman problem
% c --- 0 --- Correlation parameter

%------------------------------- Reference --------------------------------
% D. Corne and J. Knowles. Techniques for highly multiobjective
% optimisation: some nondominated points are better than others.
% Proceedings of the Annual Conference on Genetic and Evolutionary
% Computation, 2007, 773-780.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        C;  % Adjacency matrix of each map
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            c = obj.ParameterSet(0);
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 30; end
            obj.encoding = 5 + zeros(1,obj.D);
            % Randomly generate the adjacency matrices
            file = sprintf('MOTSP-M%d-D%d-c%.4f.mat',obj.M,obj.D,c);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'C');
            else
                C = cell(1,obj.M);
                C{1} = rand(obj.D);
                for i = 2 : obj.M
                    C{i} = c*C{i-1} + (1-c)*rand(obj.D);
                end
                for i = 1 : obj.M
                    C{i} = tril(C{i},-1) + triu(C{i}',1);
                end
                save(file,'C');
            end
            obj.C = C;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [sorted,rank] = sort(PopDec,2);
            index = any(sorted~=repmat(1:size(PopDec,2),size(PopDec,1),1),2);
            PopDec(index,:) = rank(index,:);
            PopObj = zeros(size(PopDec,1),obj.M);
            for i = 1 : size(PopDec,1)
                for j = 1 : obj.M
                    for k = 1 : size(PopDec,2)-1
                        PopObj(i,j) = PopObj(i,j) + obj.C{j}(PopDec(i,k),PopDec(i,k+1));
                    end
                    PopObj(i,j) = PopObj(i,j) + obj.C{j}(PopDec(i,end),PopDec(i,1));
                end
            end
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = zeros(1,obj.M) + obj.D;
        end
    end
end