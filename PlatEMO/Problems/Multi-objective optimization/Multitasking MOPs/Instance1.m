classdef Instance1 < PROBLEM
% <multi> <real> <large/none> <multitask>
% Multitasking multi-objective problem (ZDT4-R + ZDT4-G)
% SubD --- 10,10 --- Number of decision variables of each task

%------------------------------- Reference --------------------------------
% A. Gupta, Y. Ong, L. Feng, and K. C. Tan, Multiobjective multifactorial
% optimization in evolutionary multitasking, IEEE Transactions on
% Cybernetics, 2017, 47(7): 1652-1665.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        SubM;   % Number of objectives of each task
        SubD;   % Number of decision variables of each task
        L1;   	% Low bounds of the first task
        L2;   	% Low bounds of the second task
        U1;   	% Upper bounds of the first task
        U2;   	% Upper bounds of the second task
        rotmx;  % Rotation matrix
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.SubD      = obj.ParameterSet([10,10]);
            obj.SubM      = [2 2];
            obj.M         = max(obj.SubM);
            obj.D         = max(obj.SubD) + 1;
            obj.L1        = [0,zeros(1,obj.SubD(1)-1)-5];
            obj.U1        = [1,zeros(1,obj.SubD(1)-1)+5];
            obj.L2        = [0,zeros(1,obj.SubD(2)-1)-512];
            obj.U2        = [1,zeros(1,obj.SubD(2)-1)+512];
            obj.lower     = [zeros(1,obj.D-1),1];
            obj.upper     = [ones(1,obj.D-1),length(obj.SubD)];
            obj.encoding  = [ones(1,obj.D-1),2];
            ranmx         = rand(RandStream('mlfg6331_64','Seed',1),obj.D-2);
            [obj.rotmx,~] = qr(ranmx);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = nan(size(PopDec,1),obj.M);
            for i = 1 : size(PopDec,1)
                if PopDec(i,end) == 1       % Task 1
                    x1 = obj.L1 + PopDec(i,1:obj.SubD(1)).*(obj.U1-obj.L1);
                    x1(2:end) = obj.rotmx*x1(2:end)';
                    g = 1 + 10*(obj.SubD(1)-1) + sum(x1(2:end).^2-10*cos(4*pi*x1(2:end)));
                    PopObj(i,1) = x1(1);
                    PopObj(i,2) = g*(1-sqrt(x1(1)/g));
                elseif PopDec(i,end) == 2   % Task 2
                    x2 = obj.L2 + PopDec(i,1:obj.SubD(2)).*(obj.U2-obj.L2);
                    x2(2:end) = obj.rotmx*x2(2:end)';
                    g = 2 + sum(x2(2:end).^2)/4000 - prod(cos(x2(2:end)./sqrt(2:obj.SubD(2))));
                    PopObj(i,1) = x2(1);
                    PopObj(i,2) = g*(1-sqrt(x2(1)/g));
                end
            end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,N)';
            R(:,2) = 1 - R(:,1).^0.5;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            PopDec = Population.decs;
            Label  = PopDec(:,end);
            tempStream = RandStream('mlfg6331_64','Seed',2);
            for i = 1 : max(Label)
                color = rand(tempStream,1,3);
                Draw(Population(Label==i).objs+(i-1)*0.05,'o','MarkerSize',6,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                Draw(obj.PF+(i-1)*0.05,'-','LineWidth',1,'Color',color);
            end
        end
    end
end