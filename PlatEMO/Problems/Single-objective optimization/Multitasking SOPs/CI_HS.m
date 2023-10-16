classdef CI_HS < PROBLEM
% <single> <real> <large/none> <multitask>
% Multitasking problem (Griewank function + Rastrigin function)
% SubD --- 50,50 --- Number of decision variables of each task

%------------------------------- Reference --------------------------------
% K. K. Bali, Y. Ong, A. Gupta, and P. S. Tan, Multifactorial evolutionary
% algorithm with online transfer parameter estimation: MFEA-II, IEEE
% Transactions on Evolutionary Computation, 2020, 24(1): 69-83.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    properties
        SubD;   % Number of decision variables of each task
        L1;   	% Low bounds of the first task
        L2;   	% Low bounds of the second task
        U1;   	% Upper bounds of the first task
        U2;   	% Upper bounds of the second task
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.SubD     = obj.ParameterSet([50,50]);
            obj.M        = 1;
            obj.D        = max(obj.SubD) + 1;
            obj.L1       = zeros(1,obj.SubD(1)) - 600;
            obj.U1       = zeros(1,obj.SubD(1)) + 600;
            obj.L2       = zeros(1,obj.SubD(2)) - 5.12;
            obj.U2       = zeros(1,obj.SubD(2)) + 5.12;
            obj.lower    = [zeros(1,obj.D-1),1];
            obj.upper    = [ones(1,obj.D-1),length(obj.SubD)];
            obj.encoding = [ones(1,obj.D-1),2];
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),1);
            for i = 1 : size(PopDec,1)
                if PopDec(i,end) == 1       % Task 1
                    x1 = obj.L1 + PopDec(i,1:obj.SubD(1)).*(obj.U1-obj.L1);
                    PopObj(i) = 1/4000*sum(x1.^2,2) - prod(cos(x1./sqrt(repmat(1:size(x1,2),size(x1,1),1))),2) + 1;
                elseif PopDec(i,end) == 2   % Task 2
                    x2 = obj.L2 + PopDec(i,1:obj.SubD(2)).*(obj.U2-obj.L2);
                    PopObj(i) = sum(x2.^2-10*cos(2*pi*x2)+10,2);
                end
            end
        end
    end
end