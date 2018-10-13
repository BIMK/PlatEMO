function varargout = MaF2(Operation,Global,input)
% <problem> <MaF>
% A benchmark test suite for evolutionary many-objective optimization
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % This problem is DTLZ2BZ
    switch Operation
        case 'init'
            Global.M          = 3;
            Global.D          = Global.M + 9;
            Global.lower      = zeros(1,Global.D);
            Global.upper      = ones(1,Global.D);
            Global.operator   = @EAreal;
            Global.evaluation = max(1e5,1e4*Global.D);

            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,D]  = size(PopDec);
            M      = Global.M;
            
            g = zeros(N,M);
            for m = 1 : M
                if m < M
                    g(:,m) = sum(((PopDec(:,M+(m-1)*floor((D-M+1)/M):M+m*floor((D-M+1)/M)-1)/2+1/4)-0.5).^2,2);
                else
                    g(:,m) = sum(((PopDec(:,M+(M-1)*floor((D-M+1)/M):D)/2+1/4)-0.5).^2,2);
                end
            end
            PopObj = (1+g).*fliplr(cumprod([ones(size(g,1),1),cos((PopDec(:,1:M-1)/2+1/4)*pi/2)],2)).*[ones(size(g,1),1),sin((PopDec(:,M-1:-1:1)/2+1/4)*pi/2)];
            
            PopCon = [];
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            M = Global.M;
            f = UniformPoint(input,M);
            c = zeros(size(f,1),M-1);
            for i = 1 : size(f,1)
                for j = 2 : M
                    temp = f(i,j)/f(i,1)*prod(c(i,M-j+2:M-1));
                    c(i,M-j+1) = sqrt(1/(1+temp^2));
                end
            end
            if M > 5
                c = c.*(cos(pi/8)-cos(3*pi/8)) + cos(3*pi/8);
            else
                c(any(c<cos(3*pi/8)|c>cos(pi/8),2),:) = [];
            end
            f = fliplr(cumprod([ones(size(c,1),1),c(:,1:M-1)],2)).*[ones(size(c,1),1),sqrt(1-c(:,M-1:-1:1).^2)];
            varargout = {f};
    end
end