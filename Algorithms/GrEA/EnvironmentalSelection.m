function Population = EnvironmentalSelection(Population,N,div)
% The environmental selection of GrEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;

    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    Choose = LastSelection(Population(Last).objs,N-sum(Next),div);
    Next(Last(Choose)) = true;
    % Population for next generation
    Population = Population(Next);
end

function Choose = LastSelection(PopObj,K,div)
% Select part of the solutions in one front by the grid

    [N,M] = size(PopObj);

    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);
    fmin = min(PopObj,[],1);
    lb   = fmin-(fmax-fmin)/2/div;
    ub   = fmax+(fmax-fmin)/2/div;
    d    = (ub-lb)/div;
    lb   = repmat(lb,N,1);
    d    = repmat(d,N,1);
    GLoc = floor((PopObj-lb)./d);
    GLoc(isnan(GLoc)) = 0;
    
    %% Calculate GR, GCD, GCPD and GD values of each solution
    GR   = sum(GLoc,2);
    GCD  = zeros(1,N);
    GCPD = sqrt(sum(((PopObj-(lb+GLoc.*d))./d).^2,2));
    GD   = inf(N);
    for i = 1 : N-1
        for j = i+1 : N
            GD(i,j) = sum(abs(GLoc(i,:)-GLoc(j,:)));
            GD(j,i) = GD(i,j);
        end
    end

    %% Detect the grid-based dominance relation of each two solutions
    G = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(GLoc(i,:)<GLoc(j,:))-any(GLoc(i,:)>GLoc(j,:));
            if k == 1
                G(i,j) = true;
            elseif k == -1
                G(j,i) = true;
            end
        end
    end

    %% Environmental selection
    Remain = true(1,N);
    while sum(Remain) > N-K      
        % Choose the best one among the remaining solutions in the front
        CanBeChoose = find(Remain);
        temp  = find(GR(CanBeChoose)==min(GR(CanBeChoose)));
        temp2 = find(GCD(CanBeChoose(temp))==min(GCD(CanBeChoose(temp))));
        [~,q] = min(GCPD(CanBeChoose(temp(temp2))));
        q     = CanBeChoose(temp(temp2(q)));
        Remain(q)   = false;      
        % Update the GCD values
        GCD = GCD+max(M-GD(q,:),0);
        % Update the GR values
        Eq  = GD(q,:)==0 & Remain;
        Gq  = G(q,:) & Remain;
        NGq = Remain.*(1-Gq);
        Nq  = GD(q,:)<M & Remain;
        GR(Eq) = GR(Eq) + M+2;
        GR(Gq) = GR(Gq) + M;
        PD = zeros(N,1);
        for p = find((Nq.*NGq).*(1-Eq))
            if PD(p) < M-GD(q,p)
                PD(p) = M - GD(q,p);
                Gp = G(p,:) & Remain;
                for r = find(Gp.*(1-(Gq+Eq)))
                    if PD(r) < PD(p)
                        PD(r) = PD(p);
                    end
                end
            end
        end
        pp = logical(NGq.*(1-Eq));       
        GR(pp) = GR(pp) + PD(pp);
    end
    Choose = ~Remain;
end