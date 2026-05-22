classdef GDVTSF_Utils
% GDVTSF_Utils - Static class for GDVTSF specific utility and auxiliary functions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    methods(Static)
        %% 1. Calculate Fitness
        function fitness = cal_fitness(PopObj)
            N      = size(PopObj,1);
            fmax   = max(PopObj,[],1);
            fmin   = min(PopObj,[],1);
            PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
            Dis    = inf(N);
            for i = 1 : N
                SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
                for j = [1:i-1,i+1:N]
                    Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
                end
            end
            fitness = min(Dis,[],2);
        end
        
        %% 2. Calculate GDV
        function [center,idx,rd,GDV] = cal_GDV(N, D, K, T_clu, popDec, popObj, Lrd, LpopObj, Lcenter, Lidx)
            GDV = zeros(N,D);
            rd  = zeros(K,D);
            R   = rand(N,1);
            
            [idx,center] = kmeans(popDec, K, "MaxIter", T_clu);
            
            for i = 1 : K
                dis      = dist(center(i,:),Lcenter'); 
                [~,mi]   = min(dis,[],2);
                ch       = idx==i;
                C1popObj = popObj(ch,:);
                C1num    = sum(ch);
                C2popObj = LpopObj(Lidx==mi,:);
                
                fitness = GDVTSF_Utils.cal_fitness([C1popObj;C2popObj]);
                fmax    = max(fitness,[],1)+0.001;
                fmin    = min(fitness,[],1);
                fitness = (fitness-repmat(fmin,size(fitness,1),1))./repmat(fmax-fmin,size(fitness,1),1);
                
                f1      = max(fitness(1:C1num));
                f2      = max(fitness(C1num+1:end));
                rd(i,:) = mean(abs(center(i,:)-popDec(ch,:)), 1);
                
                if xor(f1 > f2, rd(i,:) < Lrd(mi,:))
                    om = 0;
                else
                    om = 1;
                end
                vec       = sign(f1-f2)*(center(i,:)-Lcenter(mi,:))+om*(center(i,:)-popDec(ch,:));
                GDV(ch,:) = R(ch,:).*rd(i,:).*vec./norm(vec);
            end
            
            if isnan(GDV)
                GDV = zeros(K,D);
            end
        end
        
        %% 3. Operator FDV
        function Offspring = operator_FDV(Problem,Rate,Acc,OffDec,OffVel)
            % Fuzzy Evolution Sub-stages Division
            Total = 1;
            S     = floor(sqrt(2*Rate*Total/Acc));
            Step  = zeros(1,S+2);  
            for i = 1 : S
                Step(i+1) = (S*i-i*i/2)*Acc;
            end
            Step(S+2) = Rate*Total;  
            
            % Fuzzy Operation
            R    = Problem.upper-Problem.lower;
            iter = Problem.FE/Problem.maxFE;  
            for i = 1 : S+1
                if iter>Step(i) && iter<=Step(i+1)
                    gamma_a = R*10^-i.*floor(10^i*R.^-1.*(OffDec-Problem.lower)) + Problem.lower;
                    gamma_b = R*10^-i.*ceil(10^i*R.^-1.*(OffDec-Problem.lower)) + Problem.lower;
                    miu1    = 1./(OffDec-gamma_a);
                    miu2    = 1./(gamma_b-OffDec);
                    logical = miu1-miu2>0;
                    OffDec  = gamma_b;
                    OffDec(logical) = gamma_a(logical);
                end
            end
            
            if nargin > 4
                Offspring = Problem.Evaluation(OffDec,OffVel);
            else
                Offspring = Problem.Evaluation(OffDec);
            end
        end
        
        %% 4. Operator GDV TSO
        function popDec = operator_GDV_TSO(Problem, population, gBest, GDV, p1, p2, p3, B, G)
            N = Problem.N;
            D = Problem.D;
            
            fitness  = GDVTSF_Utils.cal_fitness(population.objs);
            gBestDec = gBest.decs;
            popDec   = population.decs;
            popVel   = zeros(N,D);
            
            swarmN       = floor(N/3);
            [p1, p2, p3] = GDVTSF_Utils.triple_swarm_select(fitness, p1, p2, p3);
            gBestSelDec  = zeros(swarmN,D);
            
            for i = 1 : swarmN
                V      = zeros(B, D);
                H      = zeros(B, 1);
                dis    = dist(popDec(p1(i),:),gBestDec');
                [~,mi] = min(dis,[],2);
                
                % Use Global control factor G from hyperparameters
                if rand() <= G
                    gBestSelDec(i,:) = gBestDec(mi,:);
                else
                    gBestSelDec(i,:) = gBestDec(randi(length(gBest), 1),:);
                end
                
                for j = 1 : B
                    t      = rand();
                    r      = rand();
                    V(j,:) = popDec(p1(i),:) + r.*(t.*(popDec(p2(i),:)-popDec(p1(i),:)+(1-t).*(popDec(p3(i),:)-popDec(p1(i),:))));
                    L1     = norm(popDec(p1(i),:)-popDec(p2(i),:));
                    L2     = norm(popDec(p2(i),:)-popDec(p3(i),:));
                    L3     = norm(popDec(p3(i),:)-popDec(p1(i),:));
                    L      = (L1+L2+L3)/3;
                    H(j)   = log(2*pi*exp(1)*(1-exp(-norm(V(j,:)-popDec(p1(i),:)).^2/L^2)))/fitness(p1(i));
                    H(j)   = H(j) + log(2*pi*exp(1)*(1-exp(-norm(V(j,:)-popDec(p2(i),:)).^2/L^2)))/fitness(p2(i));
                    H(j)   = H(j) + log(2*pi*exp(1)*(1-exp(-norm(V(j,:)-popDec(p3(i),:)).^2/L^2)))/fitness(p3(i));
                end
                
                [~,idx]         = sort(H, 'descend');
                popVel(p2(i),:) = V(idx(1),:) - popDec(p2(i),:);
                popVel(p3(i),:) = V(idx(2),:) - popDec(p3(i),:);
            end
            
            rate = Problem.FE/Problem.maxFE;
            C1   = (1.5-rate)*repmat(rand(N,1),1,D);
            C2   = (1.5-rate)*repmat(rand(N,1),1,D);
            
            popVel(p1,:) = (gBestSelDec - popDec(p1,:));
            popVel(p2,:) = C1(p2,:).*popVel(p2,:) + C2(p2,:).*(popDec(p1,:)-popDec(p2,:));
            popVel(p3,:) = C1(p3,:).*popVel(p3,:) + C2(p3,:).*(popDec(p1,:)-popDec(p3,:));
            
            velMax = repmat((Problem.upper-Problem.lower)/1.001, N, 1);
            popVel = max(min(popVel,velMax),-velMax);
            
            try
                popDec = popDec + popVel + GDV;
            catch
            end
            
            Lower  = repmat(Problem.lower,N,1);
            Upper  = repmat(Problem.upper,N,1);
            popDec = max(min(popDec,Upper),Lower);
            
            disM  = 20;
            Site1 = repmat(rand(N,1)<0.999,1,D);
            Site2 = rand(N,D) < 1/D;
            mu    = rand(N,D);
            
            temp  = Site1 & Site2 & mu<=0.5;
            popDec(temp) = popDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                           (1-(popDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
            
            temp  = Site1 & Site2 & mu>0.5; 
            popDec(temp) = popDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                           (1-(Upper(temp)-popDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
        end
        
        %% 5. Triple Swarm Select
        function [p1,p2,p3] = triple_swarm_select(fitness, p1, p2, p3)
            Change1     = fitness(p3) > fitness(p1);
            Temp        = p1(Change1);
            p1(Change1) = p3(Change1);
            p3(Change1) = Temp;
            
            Change2     = (fitness(p2) > fitness(p1)) & (fitness(p2) > fitness(p3));
            Temp        = p1(Change2);
            p1(Change2) = p2(Change2);
            p2(Change2) = Temp;
        end
        
        %% 6. Update GBest
        function gBest = update_gbest(population, N, V, theta)
            population = population(NDSort(population.objs,1)==1);
            popObj     = population.objs;
            [NN,M]     = size(popObj);
            NV         = size(V,1);
            
            popObj = popObj - repmat(min(popObj,[],1),NN,1);
            CV     = sum(max(0,population.cons),2);
            
            cosine = 1 - pdist2(V,V,'cosine');
            cosine(logical(eye(length(cosine)))) = 0;
            gamma  = min(acos(cosine),[],2);
            
            Angle         = acos(1-pdist2(popObj,V,'cosine'));
            [~,associate] = min(Angle,[],2);
            
            Next = zeros(1,NV);
            for i = unique(associate)'
                current1 = find(associate==i & CV==0);
                current2 = find(associate==i & CV~=0);
                if ~isempty(current1)
                    APD      = (1+M*theta*Angle(current1,i)/gamma(i)).*sqrt(sum(popObj(current1,:).^2,2));
                    [~,best] = min(APD);
                    Next(i)  = current1(best);
                elseif ~isempty(current2)
                    [~,best] = min(CV(current2));
                    Next(i)  = current2(best);
                end
            end
            gBest = population(Next(Next~=0));
        end
        
        %% 7. WOF Optimizer GDVTSO
        function population = WOF_optimizer_GDV_TSO(Problem, population, t, K, B, G, T_clu)
            [V, Problem.N] = UniformPoint(Problem.N, Problem.M);    
            N              = Problem.N;
            D              = Problem.D;
            swarmN         = floor(N/3);
            gBest          = GDVTSF_Utils.update_gbest(population, swarmN, V, 0);
            
            if length(population) >= 3
                Rank = randperm(length(population),swarmN*3);
            else
                Rank = [1,1,1];
            end
            
            p1 = Rank(1:swarmN);
            p2 = Rank(swarmN+1:2*swarmN);
            p3 = Rank(2*swarmN+1:end);
            
            LpopObj         = population.objs;
            LpopDec         = population.decs;
            [Lidx, Lcenter] = kmeans(LpopDec, K, "MaxIter", T_clu);
            Lrd             = zeros(K,D);
            
            for i = 1 : K
                ch       = Lidx==i;
                Lrd(i,:) = mean(abs(Lcenter(i,:)-LpopDec(ch,:)),1);
            end
            GDV = zeros(N,D);
            cnt = 0;
            
            while cnt <= t
                cnt        = cnt + 100;
                popDec     = GDVTSF_Utils.operator_GDV_TSO(Problem, population, gBest, GDV, p1, p2, p3, B, G);
                population = Problem.Evaluation(popDec);
                popObj     = population.objs;
                [Lcenter, Lidx, Lrd, GDV] = GDVTSF_Utils.cal_GDV(N, D, K, T_clu, popDec, popObj, Lrd, LpopObj, Lcenter, Lidx);
                LpopObj                   = popObj;
            end
        end
    end
end