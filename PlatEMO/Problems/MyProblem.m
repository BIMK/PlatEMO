classdef MyProblem < PROBLEM
% 自定义多目标优化问题

    methods
        function obj = MyProblem()
            % 初始化
            obj.Global.M = 2;      % 目标数量
            obj.Global.D = 2;      % 决策变量维度
            obj.Global.lower    = [0, 0];      % 变量下界
            obj.Global.upper    = [200, 150];  % 变量上界
            obj.Global.encoding = 'real';       % 实数编码
        end
        
        function PopObj = CalObj(obj, PopDec)
            % 目标函数计算
            f1 = -(4*PopDec(:,1) + 5*PopDec(:,2));  
            f2 = -PopDec(:,1);                     
            PopObj = [f1, f2];
        end
        
        function PopCon = CalCon(obj, PopDec)
            % 约束条件计算
            g1 = -(200 - PopDec(:,1) - PopDec(:,2));
            g2 = -(200 - 1.25*PopDec(:,1) - 0.75*PopDec(:,2));
            g3 = -(150 - PopDec(:,2));
            g4 = -PopDec(:,1); 
            g5 = -PopDec(:,2); 
            
            PopCon = [g1, g2, g3, g4, g5];
        end
    end
end
