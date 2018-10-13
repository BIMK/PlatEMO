function varargout = DTLZ8(Operation,Global,input)
% <problem> <DTLZ>
% Scalable Test Problems for Evolutionary Multi-Objective Optimization
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    switch Operation
        case 'init'
            Global.M        = 3;
            Global.D        = 10*Global.M;
            Global.D        = ceil(Global.D/Global.M)*Global.M;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,D]  = size(PopDec);
            M      = Global.M;
            
            PopObj = zeros(N,M);
            for m = 1 : M
                PopObj(:,m) = mean(PopDec(:,(m-1)*D/M+1:m*D/M),2);
            end
            
            PopCon = zeros(N,M);
            PopCon(:,1:M-1) = 1 - repmat(PopObj(:,M),1,M-1) - 4*PopObj(:,1:M-1);
            if M == 2
                PopCon(:,M) = 0;
            else
                minValue    = sort(PopObj(:,1:M-1),2);
                PopCon(:,M) = 1 - 2*PopObj(:,M) - sum(minValue(:,1:2),2);
            end
            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            M = Global.M;
            if M == 2
                temp = (0:1/(input-1):1)';
                f    = [(1-temp)/4,temp];
            else
                temp = UniformPoint(input/(M-1),3);
                temp(:,3) = temp(:,3) / 2;
                temp = temp(temp(:,1)>=(1-temp(:,3))/4 & temp(:,1)<=temp(:,2) & temp(:,3)<=1/3,:);
                f    = [repmat(temp(:,2),M-1,M-1),repmat(temp(:,3),M-1,1)];
                for i = 1 : M-1
                    f((i-1)*size(temp,1)+1:i*size(temp,1),i) = temp(:,1);
                end
                gap  = sort(unique(f(:,M)));
                gap  = gap(2) - gap(1);
                temp = (1/3:gap:1)';
                f    = [f;repmat((1-temp)/4,1,M-1),temp];
                f    = unique(f,'rows');
            end
            varargout = {f};
    end
end