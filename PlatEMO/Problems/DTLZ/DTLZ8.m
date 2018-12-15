classdef DTLZ8 < PROBLEM
% <problem> <DTLZ>
% Benchmark MOP proposed by Deb, Thiele, Laumanns, and Zitzler

%------------------------------- Reference --------------------------------
% K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, Scalable test problems
% for evolutionary multiobjective optimization, Evolutionary multiobjective
% Optimization. Theoretical Advances and Applications, 2005, 105-145.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function obj = DTLZ8()
            if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 10*obj.Global.M;
            end
            obj.Global.D        = ceil(obj.Global.D/obj.Global.M)*obj.Global.M;
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            M      = obj.Global.M;
            PopObj = zeros(N,M);
            for m = 1 : M
                PopObj(:,m) = mean(PopDec(:,(m-1)*D/M+1:m*D/M),2);
            end
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            M      = obj.Global.M;
            PopObj = obj.CalObj(PopDec);
            PopCon = zeros(size(PopObj,1),M);
            PopCon(:,1:M-1) = 1 - repmat(PopObj(:,M),1,M-1) - 4*PopObj(:,1:M-1);
            if M == 2
                PopCon(:,M) = 0;
            else
                minValue    = sort(PopObj(:,1:M-1),2);
                PopCon(:,M) = 1 - 2*PopObj(:,M) - sum(minValue(:,1:2),2);
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            M = obj.Global.M;
            if M == 2
                temp = (0:1/(N-1):1)';
                P    = [(1-temp)/4,temp];
            else
                temp = UniformPoint(N/(M-1),3);
                temp(:,3) = temp(:,3) / 2;
                temp = temp(temp(:,1)>=(1-temp(:,3))/4 & temp(:,1)<=temp(:,2) & temp(:,3)<=1/3,:);
                P    = [repmat(temp(:,2),M-1,M-1),repmat(temp(:,3),M-1,1)];
                for i = 1 : M-1
                    P((i-1)*size(temp,1)+1:i*size(temp,1),i) = temp(:,1);
                end
                gap  = sort(unique(P(:,M)));
                gap  = gap(2) - gap(1);
                temp = (1/3:gap:1)';
                P    = [P;repmat((1-temp)/4,1,M-1),temp];
                P    = unique(P,'rows');
            end
        end
    end
end