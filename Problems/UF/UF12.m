classdef UF12 < PROBLEM
% <problem> <UF>
% Extended rotated DTLZ3

%------------------------------- Reference --------------------------------
% Q. Zhang, A. Zhou, S. Zhao, P. N. Suganthan, W. Liu, and S. Tiwari,
% Multiobjective optimization test instances for the CEC 2009 special
% session and competition, School of CS & EE, University of Essex, Working
% Report CES-487, 2009.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Bound;
        Lamda;
        M;
    end
    methods
        %% Initialization
        function obj = UF12()
            load('UF12_parameter','Bound','Lamda','M');
            obj.Bound = Bound;
            obj.Lamda = Lamda;
            obj.M     = M;
            obj.Global.M        = 5;
            obj.Global.D        = 30;
            obj.Global.lower    = [0,zeros(1,obj.Global.D-1)-1];
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            lamda = repmat(obj.Lamda,size(PopDec,1),1);
            z     = PopDec*obj.M';
            p     = zeros(size(z));
            temp1 = z < 0;
            temp2 = z > 1;
            p(temp1) = -z(temp1);
            p(temp2) = z(temp2) - 1;
            z(temp1) = -lamda(temp1).*z(temp1);
            z(temp2) = 1 - lamda(temp2).*(z(temp2)-1);
            psum = zeros(size(PopDec,1),5);
            for i = 1 : 5
                psum(:,i) = sqrt(sum(p(:,[1:min(6-i,4),5:end]).^2,2));
            end
            g      = 100*(26+sum((z(:,5:end)-0.5).^2-cos(20.*pi.*(z(:,5:end)-0.5)),2));
            PopObj = 2./(1+exp(-psum)).*(1+repmat(1+g,1,5).*fliplr(cumprod([ones(size(g,1),1),cos(z(:,1:5-1)*pi/2)],2)).*[ones(size(g,1),1),sin(z(:,5-1:-1:1)*pi/2)]);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            P = UniformPoint(N,5);
            P = P./repmat(sqrt(sum(P.^2,2)),1,5) + 1;
        end
    end
end