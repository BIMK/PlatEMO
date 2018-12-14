classdef MOTSP < PROBLEM
% <problem> <Combinatorial MOP>
% The multi-objective traveling salesman problem
% c --- 0 --- Correlation parameter

%------------------------------- Reference --------------------------------
% D. Corne and J. Knowles, Techniques for highly multiobjective
% optimisation: some nondominated points are better than others,
% Proceedings of the 9th Annual Conference on Genetic and Evolutionary
% Computation, 2007, 773-780.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
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
        %% Initialization
        function obj = MOTSP()
            % Parameter setting
            c = obj.Global.ParameterSet(0);
            if isempty(obj.Global.M)
                obj.Global.M = 2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 30;
            end
            obj.Global.encoding = 'permutation';
            % Randomly generate the adjacency matrices
            file = sprintf('MOTSP-M%d-D%d-c%.4f.mat',obj.Global.M,obj.Global.D,c);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'C');
            else
                C = cell(1,obj.Global.M);
                C{1} = rand(obj.Global.D);
                for i = 2 : obj.Global.M
                    C{i} = c*C{i-1} + (1-c)*rand(obj.Global.D);
                end
                for i = 1 : obj.Global.M
                    C{i} = tril(C{i},-1) + triu(C{i}',1);
                end
                save(file,'C');
            end
            obj.C = C;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            M      = obj.Global.M;
            PopObj = zeros(N,M);
            for i = 1 : M
                for j = 1 : N
                    for k = 1 : D-1
                        PopObj(j,i) = PopObj(j,i) + obj.C{i}(PopDec(j,k),PopDec(j,k+1));
                    end
                    PopObj(j,i) = PopObj(j,i) + obj.C{i}(PopDec(j,D),PopDec(j,1));
                end
            end
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
            P = zeros(1,obj.Global.M) + obj.Global.D;
        end
    end
end