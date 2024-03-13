classdef UF9 < PROBLEM
% <multi> <real> <large/none>
% Unconstrained benchmark MOP

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, S. Zhao, P. N. Suganthan, W. Liu, and S. Tiwari,
% Multiobjective optimization test instances for the CEC 2009 special
% session and competition, School of CS & EE, University of Essex, Working
% Report CES-487, 2009.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 3;
            if isempty(obj.D); obj.D = 30; end
            obj.lower    = [0,0,zeros(1,obj.D-2)-2];
            obj.upper    = [1,1,zeros(1,obj.D-2)+2];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            J1 = 4 : 3 : obj.D;
            J2 = 5 : 3 : obj.D;
            J3 = 3 : 3 : obj.D;
            Y  = X - 2*repmat(X(:,2),1,obj.D).*sin(2*pi*repmat(X(:,1),1,obj.D)+repmat(1:obj.D,size(X,1),1)*pi/obj.D);
            PopObj(:,1) = 0.5*(max(0,1.1*(1-4*(2*X(:,1)-1).^2))+2*X(:,1)).*X(:,2)   + 2*mean(Y(:,J1).^2,2);
            PopObj(:,2) = 0.5*(max(0,1.1*(1-4*(2*X(:,1)-1).^2))-2*X(:,1)+2).*X(:,2) + 2*mean(Y(:,J2).^2,2);
            PopObj(:,3) = 1-X(:,2)                                                  + 2*mean(Y(:,J3).^2,2);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,3);
            R(R(:,1)>(1-R(:,3))/4 & R(:,1)<(1-R(:,3))*3/4,:) = [];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            a = linspace(0,1,20)';
            x = a*a';
            y = a*(1-a');
            z = (1-a)*ones(size(a'));
            z(x>(1-z)/4 & x<(1-z)*3/4) = nan;
            R = {x,y,z};
        end
    end
end