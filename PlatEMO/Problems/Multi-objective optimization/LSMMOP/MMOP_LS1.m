classdef MMOP_LS1 < PROBLEM
% <2024> <multi> <real> <large/none> <multitask> <sparse/none>
% Multitasking multi-objective optimization problem (SMOP1 to SMOP8)
% SubD     --- 1000,1000 --- Number of decision variables of each task
% SubTheta ---   0.1,0.1 --- Number of decision variables of each task

%------------------------------- Reference --------------------------------
% C. Wu, Y. Tian, L. Zhang, X. Xiang, and X. Zhang, A sparsity knowledge
% transfer-based evolutionary algorithm for large-scale multitasking multi-
% objective optimization, IEEE Transactions on Evolutionary Computation,
% 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    properties
        SubM;       % Number of objectives of each task
        SubD;       % Number of decision variables of each task
        SubTheta;   % Sparsity of the Pareto set        
        L;          % Low bounds of task        
        U;          % Upper bounds of task
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.SubD,obj.SubTheta] = obj.ParameterSet([1000,1000],[0.1,0.1]);
            obj.SubM = [2 2];
            obj.M    = max(obj.SubM);
            obj.D    = max(obj.SubD) + 1;
            for i = 1:length(obj.SubD)
                obj.L{i} = [zeros(1,obj.M-1)+0,zeros(1,obj.SubD(i)-obj.M+1)-1];
                obj.U{i} = [zeros(1,obj.M-1)+1,zeros(1,obj.SubD(i)-obj.M+1)+2];
            end
            obj.lower    = [zeros(1,obj.D-1),1];
            obj.upper    = [ones(1,obj.D-1),length(obj.SubD)];
            obj.encoding = [ones(1,obj.D-1),2];
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = nan(size(PopDec,1),obj.M);            

            sk6 = find(PopDec(:,end)==1);
            X6  = PopDec(sk6,1:obj.SubD(1));
            X6  = (obj.L{1} + (obj.U{1} - obj.L{1}) .* X6) .* (X6 ~= 0);
            K6  = ceil(obj.SubTheta(1)*(obj.SubD(1)-obj.M+1));
            g66 = g4(X6(:,obj.M:end),repmat(linspace(0,1,obj.SubD(1)-obj.M+1),size(X6,1),1));
            [g66,rank] = sort(g66,2);
            temp = false(size(rank));
            for i = 1 : size(rank,1)
                temp(i,X6(i,obj.M-1+rank(i,:))==0) = true;
            end
            temp(:,1:K6) = false;
            g66(temp)    = 0;
            g66 = sum(g66,2);
            PopObj(sk6,:) = repmat(1+g66/(obj.SubD(1)-obj.M+1),1,obj.M).*fliplr(cumprod([ones(size(X6,1),1),1-cos(X6(:,1:obj.M-1)*pi/2)],2)).*[ones(size(X6,1),1),1-sin(X6(:,obj.M-1:-1:1)*pi/2)];

            sk8 = find(PopDec(:,end)==2);
            X8  = PopDec(sk8,1:obj.SubD(2));
            X8  = (obj.L{2} + (obj.U{2} - obj.L{2}) .* X8) .* (X8 ~= 0);
            K8  = ceil(obj.SubTheta(2)*(obj.SubD(2)-obj.M+1));
            g88 = sum(g3(X8(:,obj.M:obj.M+K8-1),mod(X8(:,obj.M+1:obj.M+K8)+pi,2)),2) + sum(g3(X8(:,obj.M+K8:end-1),X8(:,obj.M+K8+1:end)*0.9),2);
            PopObj(sk8,:) = repmat(1+g88/(obj.SubD(2)-obj.M+1),1,obj.M).*fliplr(cumprod([ones(size(X8,1),1),cos(X8(:,1:obj.M-1)*pi/2)],2)).*[ones(size(X8,1),1),sin(X8(:,obj.M-1:-1:1)*pi/2)];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R2 = UniformPoint(N,obj.M);
            c = ones(size(R2,1),obj.M);
            for i = 1 : size(R2,1)
                for j = 2 : obj.M
                    temp = R2(i,j)/R2(i,1)*prod(1-c(i,obj.M-j+2:obj.M-1));
                    c(i,obj.M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                end
            end
            x    = acos(c)*2/pi;
            R2   = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
            R3   = UniformPoint(N,obj.M);
            R3   = R3./repmat(sqrt(sum(R3.^2,2)),1,obj.M);      
            R{1} = [R2(1:end-1,:);1,0];
            R{2} = [R3(1:end-1,:);1,0];    
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
        %% Generate multitask indicator values
        function score = CalMetric(obj,metName,Population)
            if strcmp(metName, 'Task1_IGD') || strcmp(metName, 'Task1_HV')
                score = feval(metName,Population,obj.optimum{1});
            elseif strcmp(metName, 'Task2_IGD') || strcmp(metName, 'Task2_HV')
                score = feval(metName,Population,obj.optimum{2});
            else
                error('Incorrect metName entered');
            end
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            PopDec = Population.decs;
            Label  = PopDec(:,end);
            tempStream = RandStream('mlfg6331_64','Seed',2);
            for i = 1 : max(Label)
                color = rand(tempStream,1,3);
                Draw(obj.PF{i},'-','LineWidth',1,'Color',color);
                Draw(Population(Label==i).objs,'o','MarkerSize',6,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});                
            end
        end
    end
end

function g = g1(x,t)
    g = (x-t).^2;
end

function g = g2(x,t)
    g = 2*(x-t).^2 + sin(2*pi*(x-t)).^2;
end

function g = g3(x,t)
    g = 4-(x-t)-4./exp(100*(x-t).^2);
end

function g = g4(x,t)
    g = (x-pi/3).^2 + t.*sin(6*pi*(x-pi/3)).^2;
end