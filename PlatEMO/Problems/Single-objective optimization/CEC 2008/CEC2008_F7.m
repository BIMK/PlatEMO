classdef CEC2008_F7 < PROBLEM
% <single> <real> <large/none> <expensive/none>
% FastFractal 'DoubleDip' function

%------------------------------- Reference --------------------------------
% K. Tang, X. Yao, P. N. Suganthan, C. MacNish, Y.-P. Chen, C.-M. Chen, and
% Z. Yang, Benchmark functions for the CEC'2008 special session and
% competition on large scale global optimization, Nature Inspired
% Computation and Applications Laboratory, USTC, China, 2007.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        ran11;  % Fixed random numbers
        ran12;
        ran2;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'CEC2008.mat'),'Data');
            obj.ran11 = Data{7}(1,:);
            obj.ran12 = Data{7}(2,:);
            obj.ran2  = Data{7}(3,:);
            obj.M = 1;
            if isempty(obj.D); obj.D = 100; end
            obj.D = min(obj.D,1000);
            obj.lower    = zeros(1,obj.D) - 1;
            obj.upper    = zeros(1,obj.D) + 1;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            PopObj = zeros(size(X,1),1);
            for i = 1 : size(X,1)
                loc1 = 0;
                loc2 = 0;
                for j = 1 : size(X,2)
                    for k = 1 : 3
                        for k1 = 1 : 2^(k-1)
                            loc2 = loc2 + 1;
                            for k2 = 1 : obj.ran2(loc2)
                                loc1 = loc1 + 1;
                                if abs(X(i,j)) < 0.5
                                    PopObj(i) = PopObj(i) + (-6144*(X(i,j)-obj.ran11(loc1)).^6+3088*(X(i,j)-obj.ran11(loc1)).^4-392*(X(i,j)-obj.ran11(loc1)).^2+1)/2^(k-1)/(2-obj.ran12(loc1));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end