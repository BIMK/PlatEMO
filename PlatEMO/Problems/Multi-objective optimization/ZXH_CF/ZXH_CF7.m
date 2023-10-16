classdef ZXH_CF7 < PROBLEM
% <multi/many> <real> <large/none> <constrained>
% Constrained benchmark MOP proposed by Zhou, Xiang, and He

%------------------------------- Reference --------------------------------
% Y. Zhou, Y. Xiang, and X. He, Constrained multiobjective optimization:
% Test problem construction and performance evaluations, IEEE Transactions
% on Evolutionary Computation, 2021, 25(1): 172-186.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        k;  % Number of constrained variables
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+10;  end
            obj.lower    = zeros(1,obj.D) + 1e-10;
            obj.upper    = ones(1,obj.D)  - 1e-10;
            obj.encoding = ones(1,obj.D);
            if obj.M <= 3
                obj.k = obj.M - 1;
            elseif obj.M > 3 && obj.M <= 8 
                obj.k = floor(obj.M/2); 
            else
                obj.k = 3; 
            end
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            PopDec = varargin{1};
            OptX   = 0.2;               
            [N,D]  = size(PopDec);
            M      = obj.M;
            % Step 1: Compute cumsum 
            Sx = cumsum(PopDec(:,1:M).^2,2,'reverse'); 
            % Step 2: Compute theta
            THETA = 2/pi*atan(sqrt(Sx(:,2:end))./PopDec(:,1:M-1));
            % Step 3: Calculate Griewank function
            h = 5*(sum((PopDec(:,M+1:end)-OptX).^2,2)-prod(cos(10*pi*(PopDec(:,M+1:end)-OptX)./repmat(sqrt(1:(D-M)),N,1)),2)+1);
            % Step 4: Compute T_
            T = (1 - Sx(:,1)).^2 + h;
            % Step 5: Objectives (linear)
            G      = [ones(N,1) cumprod(THETA,2)] .* [1-THETA ones(N,1)];
            PopObj = G .* repmat((1+T),1,M) ;
            % Step 6: Constraints
            PopCon(:,1) = Sx(:,1) + h - 1; 
            PopCon(:,2) = -(Sx(:,1) + h - 1/2); 
            for i = 1 : obj.k
                PopCon(:,i+2) = max(1/4-THETA(:,i),THETA(:,i)-3/4);
            end
            Population = SOLUTION(varargin{1},PopObj,PopCon,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            T = zeros(size(R));
            for i = obj.M-1 : -1 : 1
                T(:,i) = R(:,i+1)./R(:,i)./(1-T(:,i+1)+R(:,i+1)./R(:,i));
            end
            THETA = T(:,1:obj.k);
            Valid = all(THETA>=1/4&THETA<=3/4,2);
            R     = R(Valid,:);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,1,30)';
                R = {a*a',a*(1-a'),(1-a)*ones(size(a'))};
                T2 = R{3}./R{2}./(1+R{3}./R{2});
                T1 = R{2}./(1-T2);
                THETA = cat(3,T1,T2);
                Valid = all(THETA>=1/4&THETA<=3/4,3);
                R{1}(~Valid) = nan;
            else
                R = [];
            end
        end
    end
end