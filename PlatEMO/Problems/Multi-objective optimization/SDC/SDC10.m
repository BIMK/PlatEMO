classdef SDC10 < PROBLEM
% <multi> <real> <constrained>
% Scalable high-dimensional decicsion constraint benchamrk

%------------------------------- Reference --------------------------------
% K. Qiao, J. Liang, K. Yu, C. Yue, H. Lin, D. Zhang, and B. Qu,
% Evolutionary constrained multiobjective optimization: scalable
% high-dimensional constraint benchmarks and algorithm, IEEE Transactions
% on Evolutionary Computation, 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    properties
        THETA_;  
        a; 
        CEC_Problem;        % index of high-dimension constraint function
        Distance_problem;   % index of distance function
        HCT;                % rate of unconstrained distance function variables that are used for high-dimensional constraint function
        HCT_type;           % type of transformation function
        DCT;                % rate of unconstrained distance function variables that are used for variable linkage
        DCT_type;           % type of variable linkage function
        Shape_problem;      % tyep of shape function
        b;                  % a parameter of shape function that is used to control the overlap degree between CPF and UPF
        lu;                 % upper and lower bounds of high-dimensional constraint function
        high_D_C;           % dimensions of high-dimensional constraint function
        aaa;                % used for high-dimensional constraint function
        optial_f;           % optimal objective value of high-dimensional constraints
        q1_upper;           % upper bound of shape function
        q1_lower;           % lower bound of shape function
        h;                  % objective function value of high-dimensional constraint function
        HC;                 % constraint function value of high-dimensional constraint function
        max_D_con;          % maximal number of variables that are used for transformation operator
    end
    methods
        %% Initialization
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D); obj.D = 30; end
            
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            
            obj.a = information(10);  % 10 indicates the index of problem
            obj.CEC_Problem = obj.a(1); 
            obj.Distance_problem = obj.a(2); 
            obj.HCT = obj.a(3);
            obj.HCT_type = obj.a(4); 
            obj.DCT = obj.a(5); 
            obj.DCT_type = obj.a(6); 
            obj.Shape_problem = obj.a(7); 
            obj.b = obj.a(8);

           [obj.lu,obj.high_D_C,obj.aaa,obj.optial_f]=CEC_2006_information(obj.CEC_Problem);  
           obj.max_D_con = obj.high_D_C + ceil((obj.D-obj.M-obj.high_D_C)*obj.HCT); 
           if obj.Shape_problem==1
               obj.q1_upper = [4,4];
               obj.q1_lower = [0,0];
           elseif obj.Shape_problem==2
               obj.q1_upper = [2,2];
               obj.q1_lower = [0,0];
           end
        end
        %% Calculate objective values
        function Population = Evaluation(obj,varargin)
            X  = varargin{1};
            X  = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            
            PopDec = X;
            [N,D]  = size(PopDec);
            
            % Calculate objectives and constraints of high-dimensional
            % constraint variables
            new_P = Transformation_operator(PopDec(:,obj.M+1:obj.M+obj.max_D_con),obj);
            [obj.h, obj.HC] = CEC_2006_fitness(new_P,obj.CEC_Problem,obj.aaa,obj.optial_f);
            
            % Variable linkage for distance function variables
            P = PopDec;
            control_D = ceil((obj.D-obj.high_D_C-obj.M)*obj.DCT);
            tran_D = length(D+1-control_D:D);
            if obj.DCT_type == 1
                P(:,D+1-control_D:D) = (1+repmat((1:tran_D)./length(1:tran_D),N,1)).*P(:,D+1-control_D:D)...
                    - repmat(P(:,1),1,tran_D);
            elseif obj.DCT_type == 2
                P(:,D+1-control_D:D) = (1+cos( 0.5*pi* repmat((1:tran_D)./length(1:tran_D),N,1))).*P(:,D+1-control_D:D)...
                    - repmat(P(:,1),1,tran_D);
            end
            dis_P = Distance_function(P(:,obj.M+obj.high_D_C+1:end), obj.Distance_problem);

            % Calculate angle THETA
            obj.THETA_=zeros(N,1);
            angle = atan(abs(PopDec(:,2))./PopDec(:,1));
            index = isnan(angle);
            angle(index==1) = 1;
            obj.THETA_ =2./pi.*angle;
            
            % Calculate objectives
            Pop = (obj.q1_upper - obj.q1_lower).* PopDec(:,1:2) + obj.q1_lower;
            T_ = (1-Pop(:,1).^2-Pop(:,2).^2).^2 + obj.h + dis_P;
            G_ =  [ones(N,1) cumprod(sin(pi/2*obj.THETA_),2)] .* [cos(pi/2*obj.THETA_) ones(N,1)];
            PopObj =  G_ .* repmat((1+T_),1,obj.M);

            % Calculate constraints
            if obj.Shape_problem == 1
                cc = obj.b/10;
                PopCon(:,1)=cc^2*Pop(:,1).^2+Pop(:,2).^2-cc^2;
                PopCon(:,2)=-cc*(Pop(:,1)-1)-Pop(:,2);
            else
                l      = atan(Pop(:,2)./Pop(:,1));
                index = isnan(l);
                l(index==1) = 1;
                cc = obj.b/100;
                PopCon(:,1) = Pop(:,1).^2 + Pop(:,2).^2 - (cc+0.05+0.4*sin(4*l).^16).^2;
                PopCon(:,2) = (cc-0.2*sin(4*l).^8).^2 - Pop(:,1).^2 - Pop(:,2).^2;
            end
            PopCon = [PopCon,obj.HC];
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N,~)
            R(:,1) = (0:1/(N-1):1)';
            R(:,2) = 1 - R(:,1);
            R  = R./repmat(sqrt(sum(R.^2,2)),1,2);
            cc = obj.b/100;
            invalid = (cc-0.2*sin(4*atan(R(:,2)./R(:,1))).^8).^2-R(:,1).^2-R(:,2).^2 > 0;
            while any(invalid)
                R(invalid,:) = R(invalid,:).*1.001;
                invalid = (cc-0.2*sin(4*atan(R(:,2)./R(:,1))).^8).^2-R(:,1).^2-R(:,2).^2 > 0;
            end
            theta = atan(R(:,2)./R(:,1));
            hx    = (1-sum(R.^2,2)).^2;
            R     = [cos(theta).*(1+hx),sin(theta).*(1+hx)];
            if nargin <= 2
                R = R(NDSort(R,1)==1,:);
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
        	R = obj.GetOptimum(100,[]);
            R(NDSort(R,1)>1,:) = nan;
        end
    end
end

function new_P = Transformation_operator(P,obj)
    opti_lower = obj.lu(1,:);
    opti_upper = obj.lu(2,:);
    if size(P,2) > obj.high_D_C
        chushu = floor(size(P,2)/obj.high_D_C);
        yushu = mod(size(P,2),obj.high_D_C);
        aa = cell(1,obj.high_D_C);

        for j = 1 : obj.high_D_C
            for m = 1 : chushu
                aa{j} = [aa{j}, (m-1)*obj.high_D_C+j];
            end
            if yushu >= j
                aa{j} = [aa{j}, m*obj.high_D_C+j];
            end
        end
        for j = 1 : obj.high_D_C
            if length(aa{j})>1
                if obj.HCT_type == 1
                    PP = P(:,aa{j});
                    a1 = sum(PP,2);
                    a2 = mod(a1,0.5);
                    tempa = -2*a2+1;
                    new_P(:,j) =  opti_lower(j) + (opti_upper(j) - opti_lower(j)) *  tempa;
                elseif obj.HCT_type == 2
                    PP = P(:,aa{j});
                    a1 = sum(PP,2);
                    a2 = mod(a1,0.5);
                    tempa = cos(a2*pi);
                    new_P(:,j) =  opti_lower(j) + (opti_upper(j) - opti_lower(j)) *  tempa;
                end
            else
                new_P(:,j) = opti_lower(j) +  (opti_upper(j) -  opti_lower(j)) * P(:,j);
            end
        end
    else
        P = Pop(:,2+1:2+obj.max_D_con);
        for j = 1 : obj.high_D_C
            new_P(:,j) = opti_lower(j) +  (opti_upper(j) -  opti_lower(j)) * P(:,j);
        end
    end
end