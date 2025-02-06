classdef MOEADCMT < ALGORITHM
% <2024> <multi> <real> <constrained>
% MOEA/D with competitive multitasking
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring
 
%------------------------------- Reference --------------------------------
% X. Chu, F. Ming, and W. Gong. Competitive multitasking for computational
% resource allocation in evolutionary constrained multi-objective
% optimization. IEEE Transactions on Evolutionary Computation, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [delta,nr] = Algorithm.ParameterSet(0.9,2);
            
            %% Evaluate the Population
             RMP = 0.2;      
             Nt = 2 ;      % The number of tasks
             Beta = 0.3;   % Learning phase
             Pmin = 0.1;   % Minimum selection probability
             
            %% Generate the weight vectors
             [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
             T = ceil(Problem.N/10);
            
            %% Detect the neighbours of each solution
             B = pdist2(W,W);
             [~,B] = sort(B,2);
             B = B(:,1:T);
            
            %% Initialization 
             Population{1} = Problem.Initialization();
             Population{2} = Problem.Initialization();
             Z{1}   = min(Population{1}.objs,[],1);
             Conmin = min(overall_cv(Population{1}.cons));
             Z{2}   = min(Population{2}.objs,[],1);
             A      = Population{1}; 
             rwd    =  zeros(Nt, 1);            % HR is used to store the historical rewards
             pro    =  1 / Nt * ones(1, Nt);	% Selection probability
            
            %% Optimization
             while Algorithm.NotTerminated(A) 
                CV    = sum(max(Population{1}.cons,0),2);
                fr    = length(find(CV<=0))/Problem.N;  
                sigma_obj = 0.3;
                sigma_cv = fr*10;
                Q = [];
                
                % Calculate the selection probability 
                if Problem.FE <= Beta * Problem.maxFE
                    % Stage 1: Evolution stage
                    pro  =  1 / Nt * ones(1, Nt);
                else
                    % Stage 2: Competition stage
                    if sum(rwd) ~= 0
                        pro   =   Pmin / Nt + (1 - Pmin) * rwd ./ sum(rwd);
                        pro   =   pro ./ sum(pro);
                    else
                        pro   =   1 / Nt * ones(1, Nt);
                    end
                end
                
                % Determine the a task based on the selection probability using roulette wheel method
                r = rand;
                for t = 1:Nt
                    if r <= sum(pro(1:t))
                        k = t;
                        break;
                    end
                end
                
                % For each solution
                for i = 1:Problem.N
                    % Choose the parents
                    if rand < delta
                        P = B(i,randperm(end));
                    else
                        P = randperm(Problem.N);
                    end

                   %% FCHT is the main task  
                    if k == 1   
                        % Generate an offspring
                        if rand < RMP
                            Offspring = OperatorDE(Problem,Population{k}(i),Population{k}(P(1)),Population{k}(P(2))); 
                        else
                            Offspring = OperatorDE(Problem,Population{k}(i),Population{2}(P(1)),Population{2}(P(2))); 
                        end
                        % Update the ideal point
                        Z{1} = min(Z{1},Offspring.obj);
                        Conmin = min(Conmin,overall_cv(Offspring.con));
                        Zmax  = max([Population{1}.objs;Offspring.obj],[],1);
                        Conmax = max(overall_cv([Population{1}.cons;Offspring.con]));   
                        % Update the solutions in P by Tchebycheff approach
                        g_old = max(abs(Population{1}(P).objs-repmat(Z{1},length(P),1))./W(P,:),[],2);    
                        g_new = max(repmat(abs(Offspring.obj-Z{1}),length(P),1)./W(P,:),[],2);         
                        cv_old = overall_cv(Population{1}(P).cons);   
                        cv_new = overall_cv(Offspring.con);
                        if Conmax > Conmin
                            cv_old(cv_old > 0) = (cv_old(cv_old > 0)-Conmin)/(Conmax-Conmin);
                            cv_new(cv_new > 0) = (cv_new(cv_new > 0)-Conmin)/(Conmax-Conmin);    
                        end
                        new_old = LG(g_old - g_new,[sigma_obj 0]).*max(LG(cv_old - repmat(cv_new,length(P),1),[sigma_cv 0]),0.0001);   
                        old_new = LG(g_new - g_old,[sigma_obj 0]).*max(LG(repmat(cv_new,length(P),1) - cv_old,[sigma_cv 0]),0.0001);   
                        Population{1}(P(find(new_old>=old_new,nr))) = Offspring;

                    % for task2 - CMOEAD
                        Z{2} = min(Z{2},Offspring.obj);
                        % Update the solutions in P by PBI approach
                        normW   = sqrt(sum(W(P,:).^2,2));
                        normP   = sqrt(sum((Population{2}(P).objs-repmat(Z{2},length(P),1)).^2,2));
                        normO   = sqrt(sum((Offspring.obj-Z{2}).^2,2));
                        CosineP = sum((Population{2}(P).objs-repmat(Z{2},length(P),1)).*W(P,:),2)./normW./normP;
                        CosineO = sum(repmat(Offspring.obj-Z{2},length(P),1).*W(P,:),2)./normW./normO;
                        g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                        g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                        Population{2}(P(g_old>=g_new)) = Offspring; 
                    end
                    
                   %% CMOEAD is the main task
                    if k == 2  
                        % Generate an offspring
                        if rand < RMP
                            Offspring = OperatorGAhalf(Problem,Population{k}(P(1:2))); 
                        else
                            Offspring = OperatorGAhalf(Problem,Population{1}(P(1:2))); 
                        end
                        
                        % Update the ideal point
                        Z{2} = min(Z{2},Offspring.obj);
                        % Update the solutions in P by PBI approach
                        normW   = sqrt(sum(W(P,:).^2,2));
                        normP   = sqrt(sum((Population{2}(P).objs-repmat(Z{2},length(P),1)).^2,2));
                        normO   = sqrt(sum((Offspring.obj-Z{2}).^2,2));
                        CosineP = sum((Population{2}(P).objs-repmat(Z{2},length(P),1)).*W(P,:),2)./normW./normP;
                        CosineO = sum(repmat(Offspring.obj-Z{2},length(P),1).*W(P,:),2)./normW./normO;
                        g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                        g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                        Population{2}(P(g_old>=g_new)) = Offspring; 
                        
                    % for task1 - MOEAD-FCHT
                        % Update the ideal point
                        Z{1} = min(Z{1},Offspring.obj);
                        Conmin = min(Conmin,overall_cv(Offspring.con));
                        Zmax  = max([Population{1}.objs;Offspring.obj],[],1);
                        Conmax = max(overall_cv([Population{1}.cons;Offspring.con]));
                        % Update the solutions in P by Tchebycheff approach
                        g_old = max(abs(Population{1}(P).objs-repmat(Z{1},length(P),1))./W(P,:),[],2);    
                        g_new = max(repmat(abs(Offspring.obj-Z{1}),length(P),1)./W(P,:),[],2);         
                        cv_old = overall_cv(Population{1}(P).cons);   
                        cv_new = overall_cv(Offspring.con);
                        if Conmax > Conmin
                            cv_old(cv_old > 0) = (cv_old(cv_old > 0)-Conmin)/(Conmax-Conmin);
                            cv_new(cv_new > 0) = (cv_new(cv_new > 0)-Conmin)/(Conmax-Conmin);    
                        end
                        new_old = LG(g_old - g_new,[sigma_obj 0]).*max(LG(cv_old - repmat(cv_new,length(P),1),[sigma_cv 0]),0.0001);   
                        old_new = LG(g_new - g_old,[sigma_obj 0]).*max(LG(repmat(cv_new,length(P),1) - cv_old,[sigma_cv 0]),0.0001);   
                        Population{1}(P(find(new_old>=old_new,nr))) = Offspring;
                    end
                    Q = [Q Offspring];
                end
                
                if size(Q,2) > 0
                    s = size(A,2);
                    [A,Next] =  ArchiveUpdate([A Q],Problem.N);    % update Archive
                    if s >= Problem.N     
                        rwd(k)   =  rwd(k) + sum(Next(Problem.N+1:end))/Problem.N;   % update the reward of the main task
                    end
                end 
             end
        end
    end  
end

function result = overall_cv(cv)
	cv(cv <= 0) = 0;cv = abs(cv);
	result = sum(cv,2);
end

function y = LG(x, para)
    if (x==0)
        y = repmat(0.5,length(x),1);
    else
        y = (1+exp(-1*(x-para(2))./repmat(para(1),length(x),1))).^(-1);
    end
end