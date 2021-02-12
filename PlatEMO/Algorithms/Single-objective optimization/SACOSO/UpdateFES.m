function [tArchive,SwarmFES,PbestFES,GbestFES,notEval] = UpdateFES(iter,net,SwarmFES,PbestFES,GbestFES,BU,BD,notEval,SwarmFESt1,SwarmFESt)
% Update solutions in FES-assisted swarm

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    %% Parameters in FESOperator
    c1      = 2.05;
    [c2,c3] = deal(1.025);
    PHI     = c1+c2+c3;
    KA      = 2*0.5/(2-PHI-sqrt(PHI^2-4*PHI));
    
    tArchive = [];
    [N,D]    = size(SwarmFES);
    D        = D - 1;
    %% Fitness determination
    if iter < 2
        new      = SOLUTION(SwarmFES(:,1:D));
        SwarmFES(:,end) = new.objs;
        tArchive = [tArchive,new];
    else
        % Set label
        notEST  = true(N,1);
        % Approximate the fitness of each particle
        preDecs = (SwarmFES(:,1:D)-BD)./repmat(BU-BD,N,1);
        srgtObj = sim(net,preDecs')';
        tempObj = srgtObj;
        tempFES = SwarmFES;
        
        % Calculte the distance matrix
        distMatrix = pdist2(SwarmFES(:,1:D),SwarmFES(:,1:D));
        distMatrix(logical(eye(N))) = inf;
        
        % Update fitness
        for i = 1 : N
            if notEST(i)
                SwarmFES(i,D+1) = srgtObj(i);
            end
            % Find the nearest neighbor of this particle
            [~,nbor] = min(distMatrix(i,:));
            if nbor > i
                if notEST(nbor)
                    % Fitness Estimation Strategy
%                     PositionVirtual    = SwarmFES(nbor,1:D) + KA*SwarmFESt1(nbor,1:D) + ...
%                         (1+KA*(1-c1*rand(1,D)-c2*rand(1,D)).*SwarmFESt(i,1:D))+KA*c1*rand(1,D).*PbestFES(i,1:D)+...
%                         KA*c2*rand(1,D).*GbestFES(1:D);
%                     dist = pdist2(PositionVirtual,[SwarmFES(i,1:D);SwarmFESt1(i,1:D);SwarmFESt(nbor,1:D);PbestFES(nbor,1:D);...
%                         SwarmFES(nbor,1:D);SwarmFESt1(nbor,1:D);SwarmFESt(i,1:D);PbestFES(i,1:D);GbestFES(1:D)]);
%                     Pa   = (1/dist(5)+1/dist(6)+1/dist(7)+1/dist(8)+1/dist(9))/(1/dist(1)+1/dist(2)+1/dist(3)+1/dist(4)+1/dist(9));
%                     WF   = Pa*(SwarmFES(i,D+1)/dist(1)+SwarmFESt1(i,D+1)/dist(2)+SwarmFESt(nbor,1+D)/dist(3)+PbestFES(nbor,D+1)/dist(4)+...
%                         GbestFES(D+1)/dist(9))-SwarmFESt1(nbor,1+D)/dist(6)-SwarmFESt(i,D+1)/dist(7)-PbestFES(i,D+1)/dist(8)-GbestFES(1+D)/dist(9);
%                     SwarmFES(nbor,D+1) = dist(5)*WF;
                    notEST(nbor)       = false;
                else
                    SwarmFES(nbor,D+1) = min(SwarmFES(nbor,D+1),srgtObj(nbor));
                end
            end
        end
    end
    
    %% Personal best determination for FES-assisted PSO
    if iter < 2
        update = SwarmFES(:,D+1) < PbestFES(:,D+1);
        PbestFES(update,:) = SwarmFES(update,:);
    else
        for i = 1 : N
            if SwarmFES(i,D+1) == srgtObj(i)
                if SwarmFES(i,D+1) < PbestFES(i,D+1)
                    PbestFES(i,:) = SwarmFES(i,:);
                    notEval(i)    = true;
                end
            else
                if SwarmFES(i,D+1) < PbestFES(i,D+1) && srgtObj(i) < PbestFES(i,D+1)
                    new = SOLUTION(SwarmFES(i,1:D));
                    SwarmFES(i,D+1) = new.objs;
                    if new.objs < PbestFES(i,D+1)
                        PbestFES(i,:) = SwarmFES(i,:);
                        notEval(i)    = false;
                    end
                    tArchive = [tArchive,new];
                end
            end
        end
        if isempty(tArchive)
            diff = abs(tempFES(:,1+D)-tempObj);
            DF   = mean(diff);
            select = find(diff>DF);
            if ~isempty(select)
                news = SOLUTION(SwarmFES(select,1:D));
                SwarmFES(select,1+D) = news.objs;
                tArchive = [tArchive,news];
                update   = news.objs < PbestFES(select,D+1);
                PbestFES(select(update),:) = SwarmFES(select(update),:);
                notEval(select(update))    = false;
            end  
        end
    end
    
    %% Global best position determination for FES-assisted PSO
    [value,best] = min(PbestFES(:,1+D));
    if value < GbestFES(:,1+D)
        if notEval(best)
            new  = SOLUTION(PbestFES(best,1:D));
            PbestFES(best,D+1) = new.objs;
            if new.objs < GbestFES(:,D+1)
                GbestFES = PbestFES(best,:);
            end
            tArchive = [tArchive,new];
        else
            GbestFES = PbestFES(best,:);
        end
    end
end