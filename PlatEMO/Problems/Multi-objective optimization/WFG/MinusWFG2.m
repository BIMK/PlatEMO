classdef MinusWFG2 < PROBLEM
% <2017> <multi/many> <real> <large/none> <expensive/none>
% Minus-WFG2
% K --- --- The position parameter, which should be a multiple of M-1

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

    properties(Access = private)
        K;  % Position parameter
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M + 9; end
            if rem(obj.D + obj.M, 2) == 0
                obj.D = obj.D + 1;
            end
            obj.K = obj.ParameterSet(obj.M-1);
            obj.D        = ceil((obj.D-obj.K)/2)*2 + obj.K;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = 2 : 2 : 2*obj.D;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D] = size(PopDec);
            M = obj.M;
            K = obj.K;
            L = D - K;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);

            z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);

            t1 = zeros(N,K+L);
            t1(:,1:K)     = z01(:,1:K);
            t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

            t2 = zeros(N,K+L/2);
            t2(:,1:K) = t1(:,1:K);
            % Same as <t2(:,i)=r_nonsep(t1(:,K+2*(i-K)-1:K+2*(i-K)),2)>
            t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
            % ---------------------------------------------------------

            t3 = zeros(N,M);
            for i = 1 : M-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
            end
            t3(:,M) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,M) = t3(:,M);

            h      = convex(x);
            h(:,M) = disc(x);
            PopObj = (repmat(D*x(:,M),1,M) + repmat(S,N,1).*h)*(-1);
            %PopObj = (0 + repmat(S,N,1).*h)*(-1);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            M = obj.M;
            X       = UniformPoint(N,M-1,'grid');
            X       = [X, zeros(size(X, 1), 1)];
            X(:, 1) = X(:, 1) * 0.1 + 0.9;
            R       = convex(X);
            R(:, M) = disc(X);
            R       = (repmat(2:2:2*M,size(R,1),1).*R+1)*(-1);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                a  = linspace(0,pi/2,300)';
                x  = (1-cos(a));
                y  = 1 - a*2/pi.*cos(10*a).^2;
                nd = reshape(NDSort([-x(:),-y(:)],1)==1,size(y));
                y(~nd) = nan;
                R = [-x*2-1,-y*4-1];
            elseif obj.M == 3
                a  = linspace(0,pi/2,40)';
                x  = (1-cos(a))*(1-cos(a'));
                y  = (1-cos(a))*(1-sin(a'));
                z  = 1 - a*ones(size(a'))*2/pi.*cos(10*a*ones(size(a'))).^2;
                nd = reshape(NDSort([-x(:),-y(:),-z(:)],1)==1,size(z));
                z(~nd) = nan;
                R = {-x*2-1,-y*4-1,-z*6-1};
            else
                R = [];
            end
        end
    end
end

function Output = s_linear(y,A)
    Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = r_nonsep(y,A)
    Output = zeros(size(y,1),1);
    for j = 1 : size(y,2)
        Temp = zeros(size(y,1),1);
        for k = 0 : A-2
            Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
        end
        Output = Output+y(:,j)+Temp;
    end
    Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = disc(x)
    Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end