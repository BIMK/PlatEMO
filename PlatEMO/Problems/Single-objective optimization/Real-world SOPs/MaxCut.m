classdef MaxCut < PROBLEM
% <2024> <single> <binary> <large/none>
% The max-cut problem
% dataNo --- 1 --- Number of dataset

%------------------------------- Reference --------------------------------
% Y. Tian, L. Wang, S. Yang, J. Ding, Y. Jin, and X. Zhang. Neural
% network-based dimensionality reduction for large-scale binary
% optimization with millions of variables. IEEE Transactions on
% Evolutionary Computation, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        Adj_mat;	% Adjacency matrix
        Tri_mat;	% Upper triangular matrix
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Load data
            dataNo    = obj.ParameterSet(1);
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_MC.mat'),'Dataset');
            str         = {'D941','D2344','D5000'};
            obj.Adj_mat = Dataset.(str{dataNo});
            obj.Tri_mat = sparse(triu(obj.Adj_mat, 1));
            % Parameter setting
            obj.M        = 1;
            obj.D        = length(obj.Adj_mat);
            obj.encoding = 4 + zeros(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
             PopObj = 2 * sum((sparse(PopDec) * sparse(obj.Tri_mat)) .* sparse(PopDec), 2) - PopDec * sum(obj.Tri_mat, 1)' - PopDec * sum(obj.Tri_mat, 2);
        end
    end
end