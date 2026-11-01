classdef MinusDTLZ1 < PROBLEM
% <2017> <multi/many> <real> <large/none> <expensive/none>
% Minus-DTLZ1

%------------------------------- Reference --------------------------------
% H. Ishibuchi, Y. Setoguchi, H. Masuda, and Y. Nojima. Performance of
% decomposition-based many-objective algorithms strongly depends on Pareto
% front shapes. IEEE Transactions on Evolutionary Computation, 2017, 21(2):
% 169-190.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yuma Horaguchi and Ryudai Kato (email: horaguchi-yuma-hv@ynu.jp)

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+4; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g      = 100*(obj.D-obj.M+1+sum((PopDec(:,obj.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M:end)-0.5)),2));
            PopObj = -0.5*repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(PopDec,1),1),PopDec(:,1:obj.M-1)],2)).*[ones(size(PopDec,1),1),1-PopDec(:,obj.M-1:-1:1)];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M)/2;
            g = (obj.D-obj.M+1)*220.2602;
            R = R*(1+g)*(-1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,1,10)';
                R = {a*a'/2,a*(1-a')/2,(1-a)*ones(size(a'))/2};
                g = (obj.D-obj.M+1)*220.2602;
                for i = 1:obj.M
                    R{i} = R{i}*(1+g)*(-1);
                end
            else
                R = [];
            end
        end
    end
end