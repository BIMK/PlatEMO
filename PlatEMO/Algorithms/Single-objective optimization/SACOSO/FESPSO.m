function [tArchive,SwarmFES,VelFES,PbestFES,GbestFES,notEval,notEST,PBEval] = FESPSO(iter,net,SwarmFES,VelFES,PbestFES,GbestFES,GbestRBF,Problem,notEval,notEST,PBEval,SwarmFESt1,SwarmFESt)
% FESPSO - FES-assisted particle swarm optimization

%   Example:
%       [tArchive,SwarmFES,PbestFES,GbestFES,notEval] = FESPSO(SwarmFES,VelFES,PbestFES,GbestFES,GbestRBF,Problem,notEval,SwarmFESt1,SwarmFESt)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    c1      = 2.05;
    [c2,c3] = deal(1.025);
    PHI     = 0.7298;
    [N,D]   = size(VelFES);
    BU      = Problem.upper;
    BD      = Problem.lower;
    
    %% Particle swarm optimization
    r1        = rand(N,D);
    r2        = rand(N,D);
    r3        = rand(N,D);
    VelFES    = PHI*(VelFES+c1*r1.*(PbestFES(:,1:D)-SwarmFES(:,1:D))+c2*r2.*(GbestFES(:,1:D)-SwarmFES(:,1:D))+...
        c3*r3.*(GbestRBF(:,1:D)-SwarmFES(:,1:D)));
    VelFES(:,1:D) = max(VelFES(:,1:D),repmat(BD,N,1));
    VelFES(:,1:D) = min(VelFES(:,1:D),repmat(BU,N,1));
    SwarmFES(:,1:D) = SwarmFES(:,1:D) + VelFES;
    %% Repair
    SwarmFES(:,1:D) = max(SwarmFES(:,1:D),repmat(BD,N,1));
    SwarmFES(:,1:D) = min(SwarmFES(:,1:D),repmat(BU,N,1));
    
    %% Update position
    tArchive = [];
    % Fitness determinatin
    if iter == 1
        offspring       = Problem.Evaluation(SwarmFES(:,1:D));
        SwarmFES(:,end) = offspring.objs;
        tArchive        = [tArchive,offspring];
    else
        % Set label
        FitDetermine = false(N,1);
        notEST  = true(N,1);
        % Approximate the fitness of each particle
        srgtObj = sim(net,SwarmFES(:,1:D)')';
        
        % Calculte the distance matrix
        distMatrix = pdist2(SwarmFES(:,1:D),SwarmFES(:,1:D));
        distMatrix(logical(eye(N))) = inf;
        
        % Update fitness
        [~,rankIdx] = sort(SwarmFESt(:,D+1),'ascend');
        for k = 1 : N
            i = rankIdx(k);
            FitDetermine(i) = true;
            if notEST(i)
                SwarmFES(i,D+1) = srgtObj(i);
                notEST(i)  = false;
                notEval(i) = true;
            end
            % Find the nearest neighbor of this particle
            [~,MinInd] = min(distMatrix(i,:));
            nbor       = rankIdx(MinInd);
            if notEST(nbor) || (~notEST(nbor) && notEval(nbor) && ~FitDetermine(nbor))
                % Fitness Estimation Strategy
                VirtualPosition = SwarmFES(i,1:D) + (1+PHI-PHI*c1*r1(nbor,1:D)-PHI*c2*r2(nbor,1:D)-...
                    PHI*c3*r3(nbor,1:D)).*SwarmFESt(nbor,1:D) + SwarmFESt1(i,1:D) + PHI*c1*r1(i,1:D).*PbestFES(nbor,1:D)...
                    + PHI*c2*r2(nbor,1:D).*GbestFES(1:D) + PHI*c3*r3(nbor,1:D).*GbestRBF(1:D);
                Dist = pdist2(VirtualPosition,[SwarmFES(i,1:D);SwarmFESt(nbor,1:D);SwarmFESt1(i,1:D);...
                    PbestFES(nbor,1:D);GbestFES(1:D);GbestRBF(1:D);SwarmFES(nbor,1:D);SwarmFESt(i,1:D);...
                    SwarmFESt1(nbor,1:D);PbestFES(i,1:D)]);
                Dist = [Dist,Dist(5),Dist(6)];
                if all(Dist~=0)
                    DistTemp1 = sum(Dist(1:6).^(-1));
                    DistTemp2 = sum(Dist(7:12).^(-1));
                    DistRadio = DistTemp1/DistTemp2;
                    VirtualFitness = SwarmFES(i,D+1)/Dist(1) + SwarmFESt(nbor,D+1)/Dist(2) + SwarmFESt1(i,D+1)/Dist(3) + ...
                        PbestFES(nbor,1+D)/Dist(4) + GbestFES(D+1)/Dist(5) + GbestRBF(D+1)/Dist(6);
                    if notEST(nbor)
                        SwarmFES(nbor,D+1) = Dist(7)*(VirtualFitness*DistRadio-(SwarmFESt(i,D+1)/Dist(8)+SwarmFESt1(nbor)/Dist(9)+PbestFES(i,D+1)/Dist(10)+GbestFES(D+1)/Dist(11)+GbestRBF(D+1)/Dist(12)));
                    else
                        TempFit = Dist(7)*(VirtualFitness*DistRadio-(SwarmFESt(i,D+1)/Dist(8)+SwarmFESt1(nbor,D+1)/Dist(9)+PbestFES(i,D+1)/Dist(10)+GbestFES(D+1)/Dist(11)+GbestRBF(D+1)/Dist(12)));
                        SwarmFES(nbor,D+1) = min(SwarmFES(nbor,D+1),TempFit);
                    end
                    notEST(nbor)  = false;
                    notEval(nbor) = true;
                end
            end
            Ind = find(distMatrix(i,:)==0);
            for j = 1 : Ind
                if FitDetermine(Ind(j)) == 0
                    SwarmFES(Ind(j),D+1) = SwarmFES(nbor,D+1);
                end
            end
        end
    end
    
    % Personal best determination
    if iter == 1
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
                    offspring = Problem.Evaluation(SwarmFES(i,1:D));
                    SwarmFES(i,D+1) = offspring.objs;
                    if offspring.objs < PbestFES(i,D+1)
                        PbestFES(i,:) = SwarmFES(i,:);
                        notEval(i)    = false;
                        notEST(i)     = false;
                    end
                    tArchive = [tArchive,offspring];
                end
            end
        end
        if isempty(tArchive)
            remain = find(SwarmFES(:,1+D) ~= srgtObj);
            
            for k = 1 : length(remain)
                diff = abs(SwarmFES(remain,1+D)-srgtObj(remain));
                DF   = mean(diff);
                select    = find(diff>DF);
                if ~isempty(select)
                    offspring = Problem.Evaluation(SwarmFES(remain(select),1:D));
                    SwarmFES(remain(select),1+D) = offspring.objs;
                    tArchive = [tArchive,offspring];
                    notEval(remain(select)) = false;
                    notEST(remain(select))  = false; 
                end
            end
        end
        
        % Update Personal best solution
        update = SwarmFES(:,1+D) < PbestFES(:,1+D);
        PbestFES(update,:) = SwarmFES(update,:);
        PBEval(update) = notEval(update);
    end
    
    % Global best position determination
    [value,best] = min(PbestFES(:,1+D));
    if iter == 1
        if value < GbestFES(D+1)
            GbestFES = PbestFES(best,:);
        end
    else
        if value < GbestFES(:,1+D)
            if ~PBEval(best)
                new  = Problem.Evaluation(PbestFES(best,1:D));
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
end