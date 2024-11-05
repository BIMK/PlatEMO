classdef LSCM1 < PROBLEM
% <2024> <multi> <real> <large/none> <constrained>
% Large-scale constrained multiobjective benchmark problem

%------------------------------- Reference --------------------------------
% K. Qiao, J. Liang, K. Yu, W. Guo, C. Yue, B. Qu, and P. N. Suganthan,
% Benchmark problems for large-scale constrained multi-objective
% optimization with baseline results, Swarm and Evolutionary Computation,
% 2024, 86: 101504.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao
% If you have any question, please email qiaokangjia@yeah.net

    properties
        ku_index = 1; % index of LSCM function
        base; % the dimension of each segment
        duan; % the number of segments
        DCT; % the form of unconstraint variable linkage
        CCT; % the form of constraint variable linkage
        Dis_function; % index of distance functions
        Q_form; % index of Q functions
        PF_form; % index of PF functions
        CEC_Problem; % index of constraint functions
        lu; % upper and lower bounds of constraint functions
        high_D_C; % dimensions of constraint functions
        aaa; % used for constraint functions
        optial_f; % optimal objective values of constraint functions
        
        nk; % Number of subcomponents in each subset
        sublen;	% Number of variables in each subset
        len; % Cumulative sum of lengths of variable groups
        
        CV; % store the index of constaint variables
        Dis; % store the index of distance function variables
        
        h;  % objective function value of constraint function
        HC; % constraint function value of constraint function
    end
    methods
        %% Initialization
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D)
                obj.D = 100;
                obj.base = 100;
                obj.duan =  obj.D/obj.base; % obj.duan is the parameter k in Fig. 1 of the manuscript, that is, the number of segments
            else
                obj.base = ceil(obj.D/1000) * 100;
                obj.D = round(obj.D/obj.base)*obj.base;
                obj.duan =  obj.D/obj.base;
            end
            
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            
            [a, b, obj.CCT, obj.DCT,obj.Dis_function,obj.Q_form,obj.PF_form] = duan_yueshu(obj.ku_index); % b indicates the index of CEC2006 functions
            obj.nk = a(1:obj.duan); % extract the nk values of former k segments 
            obj.CEC_Problem = b(1:obj.duan); % extract the CEC2006 constraint function index of former k segments 
            
            for j = 1 : obj.duan
                [obj.lu{j},obj.high_D_C(j),obj.aaa{j},obj.optial_f(j)] = CEC_2006_information(obj.CEC_Problem(j));
            end
            
            % Calculate the number of variables in each set
            c = 3.8*0.1*(1-0.1);
            for i = 1 : obj.M-1
                c = [c,3.8.*c(end).*(1-c(end))];
            end
            
            [obj] = variable_information(obj,c); % extract the index of differnet varaibles from the varaible vector
        end
        %%
        function Population = Evaluation(obj,varargin)
            X = varargin{1};
            X = max(min(X,repmat(obj.upper,size(X,1),1)),repmat(obj.lower,size(X,1),1));
            PopDec = X;
            [N,D]  = size(PopDec);
            M      = obj.M;
            [PopDec] = variable_linkage(PopDec,obj); 

            G = zeros(N,M);
            Con= zeros(N,1);
            for mm = 1: obj.duan
                P = PopDec(:, obj.CV{mm,1}); 
                new_P =  repmat(obj.lu{mm}(1,:),N,1) + repmat((obj.lu{mm}(2,:) - obj.lu{mm}(1,:)),N,1) .*  P;
                [h, hc] = CEC_2006_fitness(new_P, obj.CEC_Problem(mm), obj.aaa{mm},obj.optial_f(mm));
                
                for i = 1 : M
                    G(:,i) = G(:,i) + h/obj.M;
                end
                Con = Con + hc;
                
                for i = 1 : 2 : M
                    temp = 0;
                    for j = 1 : obj.nk(mm)
                        temp = temp + Distance_function( PopDec(:, obj.Dis{mm,i}{j,1}), obj.Dis_function(1));
                    end
                    G(:,i)  = G(:,i) + temp./repmat(obj.sublen{mm}(i),N,1)./obj.nk(mm);
                end
                for i = 2 : 2 : M
                    temp = 0;
                    for j = 1 : obj.nk(mm)
                        temp = temp + Distance_function( PopDec(:, obj.Dis{mm,i}{j,1}), obj.Dis_function(2));
                    end
                    G(:,i)  = G(:,i) + temp./repmat(obj.sublen{mm}(i),N,1)./obj.nk(mm);
                end
            end
            
            obj.HC = Con;
            
            PopObj = PF_function(PopDec,G,obj.PF_form,obj.Q_form,N,M);
            PopCon = Con;
            Population  = SOLUTION(X,PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            if obj.PF_form == 1
                R = UniformPoint(N,obj.M);
            elseif obj.PF_form == 2
                R = UniformPoint(N,obj.M);
                R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            elseif obj.PF_form == 3
                interval     = [0,0.251412,0.631627,0.859401];
                median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
                X            = UniformPoint(N,obj.M-1,'grid');
                X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
                X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
                R            = [X,2*(obj.M-sum(X/2.*(1+sin(3*pi.*X)),2))];
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.PF_form == 1
                if obj.M == 2
                    R = obj.GetOptimum(100);
                elseif obj.M == 3
                    a = linspace(0,1,10)';
                    R = {a*a',a*(1-a'),(1-a)*ones(size(a'))};
                else
                    R = [];
                end
            elseif obj.PF_form == 2
                if obj.M == 2
                    R = obj.GetOptimum(100);
                elseif obj.M == 3
                    a = linspace(0,pi/2,10)';
                    R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
                else
                    R = [];
                end
            elseif obj.PF_form == 3
                if obj.M == 2
                    x      = linspace(0,1,100)';
                    y      = 2*(2-x/2.*(1+sin(3*pi*x)));
                    nd     = NDSort([x,y],1)==1;
                    x(~nd) = nan;
                    R      = [x,y];
                elseif obj.M == 3
                    [x,y]  = meshgrid(linspace(0,1,20));
                    z      = 2*(3-x/2.*(1+sin(3*pi*x))-y/2.*(1+sin(3*pi*y)));
                    nd     = reshape(NDSort([x(:),y(:),z(:)],1)==1,size(z));
                    z(~nd) = nan;
                    R      = {x,y,z};
                else
                    R = [];
                end
            end
        end
    end
end