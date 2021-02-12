classdef LIRCMOP14 < PROBLEM
% <multi> <real> <large/none> <constrained>
% Constrained benchmark MOP with large infeasible regions

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Huang, Y. Fang, Y. You, J. Mo, C. Wei, and E.
% Goodman, An improved epsilon constraint-handling method in MOEA/D for
% CMOPs with large infeasible regions, Soft Computing, 2019, 23:
% 12491-12510.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenji Li
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 3;
            if isempty(obj.D); obj.D = 30; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            variable_length = size(X,2);
            popsize         = size(X,1);
            sum1            = zeros(popsize,1);
            for j = 3 : variable_length
                sum1 = sum1+10*(X(:,j)-0.5).^2;
            end
            PopObj(:,1) = (1.7057+sum1).*cos(0.5*pi*X(:,1)).*cos(0.5*pi*X(:,2));
            PopObj(:,2) = (1.7057+sum1).*cos(0.5*pi*X(:,1)).*sin(0.5*pi*X(:,2));
            PopObj(:,3) = (1.7057+sum1).*sin(0.5*pi*X(:,1));
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj      = obj.CalObj(X);
            gx          = PopObj(:,1).^2+PopObj(:,2).^2+PopObj(:,3).^2;
            PopCon(:,1) = (gx-9).*(4-gx);
            PopCon(:,2) = (gx-3.61).*(3.24-gx);
            PopCon(:,3) = (gx-3.0625).*(2.56-gx);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,3);
            R = sqrt(3.0625)*R./repmat(sqrt(sum(R.^2,2)),1,3);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            a = linspace(0,pi/2,10)';
            R = {sin(a)*cos(a')*sqrt(3.0625),sin(a)*sin(a')*sqrt(3.0625),cos(a)*ones(size(a'))*sqrt(3.0625)};
        end
    end
end