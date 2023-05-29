classdef RBFNNPC < handle
% Radial Basis Network

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        spread   = 0;
        srgtSRGT = cell(1,1);
    end
    methods
        %% Constructor
        function obj = RBFNNPC(spread)
            obj.spread = spread;
        end

        %% Train
        function train(obj,X,distanceD)
            % Input ---(x1£¬x2) £¬decision number is changed to 2*D
            % Output --1 or 0 , meaning is better or worse
            N = size(X,1);
            input = zeros(N*N,2*distanceD);
            output = zeros(N*N,1);
            FrontNo = X(:,end);
            for i = 1:N
                for j = 1:N
                    input(N*(i-1)+j,:) = [X(i,1:distanceD),X(j,1:distanceD)];
                    if FrontNo(i)<FrontNo(j)
                        output(N*(i-1)+j) = 2;
                    else
                        output(N*(i-1)+j) = 1;
                    end
                end
            end
            X = input;
            T = output;
            E = eye(N);
            Find = E(:)==1;%delete£¨x£¬x£©
            X(Find,:) = [];
            T(Find) = [];

            Tc = ind2vec(T');
            obj.srgtSRGT = newpnn(X',Tc,obj.spread);
        end

        %% Predict for data
        function Y = lastpredict(obj,X,distanceD,Preference,flag)
            % change the input including the preference point
            N = size(X,1);
            decs = X(:,1:distanceD);
            numberP = size(Preference,1);

            Pref = zeros(N,distanceD);
            for i = 1:N
                Pref(i,:) = Preference(mod(i+numberP,numberP)+1,1:distanceD);% better than a random reference point
            end

            % forward forecasting
            X  = [decs,Pref];
            Yc = sim(obj.srgtSRGT,X');
            Y = vec2ind(Yc);
            Y = Y';

            % reverse forecasting
            x=[Pref,decs];
            yc = sim(obj.srgtSRGT,x');
            y = vec2ind(yc);
            y = y';

            result = Y-y;
            if flag == 0
                Y = result;
            else
                Y(result == 0) = 1.5;
            end
        end
    end
end