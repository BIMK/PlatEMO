function varargout = WFG3(Operation,Global,input)
% <problem> <WFG>
% A Review of Multi-objective Test Problems and a Scalable Test Problem
% Toolkit
% K ---    --- The position parameter, which should be a multiple of M-1
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    K = Global.ParameterSet(Global.M-1);
    switch Operation
        case 'init'
            Global.M        = 3;
            Global.D        = Global.M + 9;
            Global.D        = ceil((Global.D-Global.M+1)/2)*2 + Global.M - 1;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = 2 : 2 : 2*Global.D;
            Global.operator = @EAreal;
            
            PopDec    = rand(input,Global.D).*repmat(2:2:2*Global.D,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,D]  = size(PopDec);
            M      = Global.M;
            
            L = D - K;
            D = 1;
            S = 2 : 2 : 2*M;
            A = [1,zeros(1,M-2)];

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

            h      = linear(x);
            PopObj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

            PopCon = [];

            varargout = {input,PopObj,PopCon};
        case 'PF'
            PopDec    = (0:1/(input-1):1)';
            PopDec    = [PopDec,zeros(input,Global.M-2)+0.5,zeros(input,1)];
            h         = linear(PopDec);
            h         = repmat(2:2:2*Global.M,size(h,1),1).*h;
            varargout = {h};
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

function Output = linear(x)
    Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end