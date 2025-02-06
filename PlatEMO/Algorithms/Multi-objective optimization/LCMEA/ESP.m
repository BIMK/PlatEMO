function Offspring = ESP(Problem, Population)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Number of promising directions
    RefNo  = 10;
    SubN   = min(ceil(Problem.N/RefNo), Problem.N);

    %% Population info
    CV   = max(Population.cons, 0);
    Cons = sum(CV, 2);
    FNum = sum(Cons==0);
    FIndex = find(Cons==0);
    
    %% Decision boundary
    Lower = repmat(Problem.lower, Problem.N, 1);
    Upper = repmat(Problem.upper, Problem.N, 1);

    %% Select starting and ending solutions
    Fitness = AdaFitness(Population.objs, Cons);
    [~, Index] = sort(Fitness, 'ascend');
    if rand() < exp(-RefNo*FNum/Problem.N) 
        NDIndex = find(Fitness<=1);
        %% Starting solutions
        if FNum < RefNo
            SIndex = Index(1:RefNo);
        else
            SIndex  = FIndex(randperm(FNum, RefNo));
        end

        %% Ending solutions
        EIndex  = Index(floor(Problem.N/2)+randperm(Problem.N-floor(Problem.N/2), RefNo));

        if FNum > 1
            RIndex = NDIndex;
        else
            RIndex = Index(1:RefNo);
        end
    else
        [FrontNo, ~] = NDSort(Population.objs, Population.cons, inf);
        NIndex = find(FrontNo==1);
        DIndex = find(FrontNo>1);
        NNum   = numel(NIndex);
        DNum   = numel(DIndex);
        if NNum < 2*RefNo || DNum < RefNo
            %% Convergence-biased selection
            if NNum < RefNo
                SIndex = Index(1:RefNo);
                EIndex = Index(floor(Problem.N/2)+randperm(Problem.N-floor(Problem.N/2), RefNo));
            else
                if DNum < RefNo
                    Select = randperm(NNum, 2*RefNo-DNum);
                    SIndex = NIndex(Select(1:RefNo));
                    EIndex = [DIndex, NIndex(Select(RefNo+1:end))];
                else
                    SIndex = NIndex(randperm(NNum, RefNo));
                    EIndex = DIndex(randperm(DNum, RefNo));
                end
            end
        else
            if rand() < exp(-0.5*RefNo*NNum/Problem.N)
                %% Diversity-biased selection
                Select = randperm(NNum, 2*RefNo);
                SIndex = NIndex(Select(1:RefNo));
                EIndex = NIndex(Select(RefNo+1:end));

            else
                %% Convergence-biased selection
                SIndex = NIndex(randperm(NNum, RefNo));
                EIndex = DIndex(randperm(DNum, RefNo));
            end 
        end
        if FNum > 1
            RIndex = NIndex;
        else
            RIndex = SIndex;
        end
    end

    PopDec = [];
    Vector = Population(EIndex).decs - Population(SIndex).decs;
    Direct = Vector./repmat(sum(Vector.^2,2).^(1/2),1,Problem.D);
    % Sampling in the promising direction
    for i = 1 : 1 : RefNo
        lambda = (Population(RIndex).decs - repmat(Population(SIndex(i)).decs,numel(RIndex),1))*Direct(i,:)';
        beta   = std(lambda);
        r = normrnd(0, beta, [SubN,1])*(1/(1+((Problem.N-FNum+1)/Problem.N))^ 2);
        OffDec = repmat(r,1,Problem.D).*repmat(Direct(i,:),SubN,1)+repmat(Population(SIndex(i)).decs,SubN,1);
        PopDec = [PopDec;OffDec];
    end
    
    % Perturbation
    BU = max(Population.decs, [], 1);
    BD = min(Population.decs, [], 1);
    delta =  ((Problem.maxFE-Problem.FE+1)/Problem.maxFE)^ 2;
    if delta < 1e-1
        delta = 0;
    end

    p = ceil(0.1*Problem.N);
    Pindex = randperm(Problem.N,Problem.N);
    
    % Scaling: perturbation % For LSCM
    if rand() < 0.1
        RNum  = Problem.N - 2*p;
        PopDec(Pindex(1+RNum:RNum+p),:) = PopDec(Pindex(1+RNum:RNum+p),:)+repmat(randn(p,1),1,Problem.D).*PopDec(Pindex(1+RNum:RNum+p),:);
        [~, Index] = sort(Fitness, 'ascend');
        PopDec2 = PopDec(Pindex(1+RNum+p:end),:);
        F = 0.5;
        ParDec = Population(randperm(Problem.N,p)).decs;
        SiteS  = rand(p,Problem.D) < repmat(rand(p,1),1,Problem.D);
        BDec   = repmat(Population(Index(randperm(RefNo,1))).decs, p, 1);
        PopDec2(SiteS) = PopDec2(SiteS) + F.*(BDec(SiteS)-ParDec(SiteS));
        PopDec(Pindex(1+RNum+p:end),:) = PopDec2;
    else
        p = 2*p;
        RNum  = Problem.N - p;
        [~, Index] = sort(Fitness, 'ascend');
        PopDec2 = PopDec(Pindex(1+RNum:end),:);
        F = 0.5;
        ParDec = Population(randperm(Problem.N,p)).decs;
        SiteS  = rand(p,Problem.D) < repmat(rand(p,1),1,Problem.D);
        BDec   = repmat(Population(Index(randperm(RefNo,1))).decs, p, 1);
        PopDec2(SiteS) = PopDec2(SiteS) + F.*(BDec(SiteS)-ParDec(SiteS));
        PopDec(Pindex(1+RNum:end),:) = PopDec2;
    end

    % Region: local perturbation
    SiteR = rand(RNum, Problem.D) < 2*repmat(rand(RNum,1), 1, Problem.D);
    PR    = repmat((BU-BD)/RefNo, RNum, 1) .* repmat(2*randn(RNum,1),1,Problem.D) * delta;
    PopDec1 = PopDec(Pindex(1:RNum),:);
    PopDec1(SiteR) = PopDec1(SiteR) + PR(SiteR);
    PopDec(Pindex(1:RNum),:) = PopDec1;

    PopDec = max(min(PopDec,Upper),Lower);
    [proM,disM] = deal(1,20);
    OffspringT = PopDec;
       
    % Polynomial mutation
    Site  = rand(Problem.N,Problem.D) < proM/Problem.D;
    mu    = rand(Problem.N,Problem.D);
    temp  = Site & mu<=0.5;
    OffspringT       = min(max(OffspringT,Lower),Upper);
    OffspringT(temp) = OffspringT(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(OffspringT(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    OffspringT(temp) = OffspringT(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-OffspringT(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring = Problem.Evaluation(OffspringT);
end