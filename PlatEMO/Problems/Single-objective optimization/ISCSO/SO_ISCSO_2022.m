classdef SO_ISCSO_2022 < PROBLEM
% <2023> <single> <integer> <large> <constrained> 
% International student competition in structural optimization

%------------------------------- Reference --------------------------------
% S. K. Azad and S. K. Azad. A standard benchmarking suite for structural
% optimization algorithms: ISCSO 2016-2022. Structures, 2023, 58, 105409.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This functions is written by Saeid Kazemzadeh Azad
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 1;
            obj.D = 336;     
            obj.lower    = ones(1,obj.D);
            obj.upper    = 37*ones(1,obj.D);
            obj.encoding = 3*ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            PopObj = zeros(size(X,1),1);
            PopCon = zeros(size(X,1),2);
            for i = 1 : size(X,1)
                [PopObj(i,1),PopCon(i,1),PopCon(i,2)] = ISCSO_2022(X(i,:),0); 
            end
            Population = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
    end
end