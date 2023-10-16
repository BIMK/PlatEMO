classdef TSP < PROBLEM
% <single> <permutation> <large/none>
% The traveling salesman problem

%------------------------------- Reference --------------------------------
% D. Corne and J. Knowles, Techniques for highly multiobjective
% optimisation: some nondominated points are better than others,
% Proceedings of the Annual Conference on Genetic and Evolutionary
% Computation, 2007, 773-780.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        R;  % Locations of points
        C;  % Adjacency matrix
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            obj.M = 1;
            if isempty(obj.D); obj.D = 30; end
            obj.encoding = 5 + zeros(1,obj.D);
            % Randomly generate the adjacency matrix
            file = sprintf('TSP-D%d.mat',obj.D);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'R');
            else
                R = rand(obj.D,2);
                save(file,'R');
            end
            obj.R = R;
            obj.C = pdist2(obj.R,obj.R);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                for j = 1 : size(PopDec,2)-1
                    PopObj(i) = PopObj(i) + obj.C(PopDec(i,j),PopDec(i,j+1));
                end
                PopObj(i) = PopObj(i) + obj.C(PopDec(i,end),PopDec(i,1));
            end
        end
        %% Display a population in the decision space
        function DrawDec(obj,Population)
            [~,best] = min(Population.objs);
            Draw(obj.R(Population(best).dec([1:end,1]),:),'-k','LineWidth',1.5);
            Draw(obj.R);
        end
    end
end