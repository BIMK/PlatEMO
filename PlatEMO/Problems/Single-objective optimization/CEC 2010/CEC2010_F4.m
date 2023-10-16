classdef CEC2010_F4 < PROBLEM
% <single> <real> <constrained>
% CEC'2010 constrained optimization benchmark problem

%------------------------------- Reference --------------------------------
% R. Mallipeddi and P. N. Suganthan, Problem definitions and evaluation
% criteria for the CEC 2010 competition on constrained real-parameter
% optimization, Nanyang Technological University, Singapore, 2010.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        O;  % Optimal decision vector
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'CEC2010.mat'),'Data');
            obj.O = Data{4};
            obj.M = 1;
            if isempty(obj.D); obj.D = 10; end
            obj.D = min(obj.D,length(obj.O));
            obj.lower    = zeros(1,obj.D) - 50;
            obj.upper    = zeros(1,obj.D) + 50;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            Z = PopDec - repmat(obj.O(1:size(PopDec,2)),size(PopDec,1),1);
            PopObj = max(Z,[],2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            Z = PopDec - repmat(obj.O(1:size(PopDec,2)),size(PopDec,1),1);
            PopCon(:,1) = abs(mean((Z.*cos(sqrt(abs(Z)))),2)) - 1e-4;
            PopCon(:,2) = abs(sum((Z(:,1:end/2-1)-Z(:,2:end/2)).^2,2)) - 1e-4;
            PopCon(:,3) = abs(sum((Z(:,end/2+1:end-1).^2-Z(:,end/2+2:end)).^2,2)) - 1e-4;
            PopCon(:,4) = abs(sum(Z,2)) - 1e-4;
        end
    end
end