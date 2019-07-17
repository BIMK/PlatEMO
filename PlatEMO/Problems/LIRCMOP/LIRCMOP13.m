classdef LIRCMOP13 < PROBLEM
% <problem> <LIR-CMOP>
% Constrained benchmark MOP with large infeasible regions

%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Huang, Y. Fang, Y. You, J. Mo, C. Wei, and E.
% Goodman, An improved epsilon constraint-handling method in MOEA/D for
% CMOPs with large infeasible regions, Soft Computing, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenji Li
    
    methods
        %% Initialization
        function obj = LIRCMOP13()
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 30;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            variable_length = size(X,2);
            popsize         = size(X,1);
            sum            = zeros(popsize,1);
            for j=3:variable_length
                sum=sum+10*(X(:,j)-0.5).^2;
            end
            PopObj(:,1)=(1.7057+sum).*cos(0.5*pi*X(:,1)).*cos(0.5*pi*X(:,2));
            PopObj(:,2)=(1.7057+sum).*cos(0.5*pi*X(:,1)).*sin(0.5*pi*X(:,2));
            PopObj(:,3)=(1.7057+sum).*sin(0.5*pi*X(:,1));
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            PopObj      = obj.CalObj(X);
            gx          =  PopObj(:,1).^2+PopObj(:,2).^2+PopObj(:,3).^2;
            PopCon(:,1) = (gx-9).*(4-gx);
            PopCon(:,2) = (gx-3.61).*(3.24-gx);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'LIRCMOP_PF.mat'),'PF');
            P = PF{13};
        end 
    end
end