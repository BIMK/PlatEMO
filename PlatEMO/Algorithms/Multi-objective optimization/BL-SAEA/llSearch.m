function [ best_xl, best_f,best_con,llFunctionEvaluations,ulFunctionEvaluations] = llSearch(xl0,flag,ulPop,P_llPop,Problem,llDim, llDimMin, llDimMax, xu,llPopSize,llStoppingCriteria,proC,disC,proM,disM,llFunctionEvaluations,llMaxGens,ulFunctionEvaluations )

    if flag==0

        %Lower population initialization (LHS sampling)
        Pl = UniformPoint(llPopSize,llDim,'Latin');
        llPop    = repmat(llDimMax-llDimMin,llPopSize,1).*Pl+repmat(llDimMin,llPopSize,1);

    else
        %random initialization
        if(llDim<=3)
            P = UniformPoint(13,llDim,'Latin');
            llPop    = repmat(llDimMax-llDimMin,13,1).*P+repmat(llDimMin,13,1);
        else
            Pl = UniformPoint(23,llDim,'Latin');
            llPop    = repmat(llDimMax-llDimMin,23,1).*Pl+repmat(llDimMin,23,1);
        end

        %Find the closest point in the upper population and the corresponding optimal solution in the lower level
        [ulPopSize,~]=size(ulPop);
        for e=1:ulPopSize
            c_o_s(e,1)=1-pdist([xu;ulPop(e,:)],'cosine');
        end
        [maxv,maxi] = max(c_o_s);
        Y=P_llPop(maxi,:);

        %disturbance
        sigma = 0.2*(llDimMax-llDimMin);
        for i=1:1
            mu       = rand(1,llDim) < 0.5;
            Ydec     = Y;
            Ydec(mu) = Ydec(mu) + sigma(mu).*randn(1,sum(mu));
            %Out of bounds processing
            for j=1:llDim
                if Ydec(j)<llDimMin(j)
                    Ydec(j)=llDimMin(j);

                elseif Ydec(j)>llDimMax(j)
                    Ydec(j)=llDimMax(j);
                end
            end

            llPop=[llPop;Ydec];
        end
        llPop=[llPop;Y];

        %disturbance
        if(llDim<=3)
            for i=1:14
                mu       = rand(1,llDim) < 0.5;
                Ydec     = xl0;
                Ydec(mu) = Ydec(mu) + sigma(mu).*randn(1,sum(mu));
                %Out of bounds processing
                for j=1:llDim
                    if Ydec(j)<llDimMin(j)
                        Ydec(j)=llDimMin(j);

                    elseif Ydec(j)>llDimMax(j)
                        Ydec(j)=llDimMax(j);
                    end
                end

                llPop=[llPop;Ydec];
            end
        else
            for i=1:24
                mu       = rand(1,llDim) < 0.5;
                Ydec     = xl0;
                Ydec(mu) = Ydec(mu) + sigma(mu).*randn(1,sum(mu));
                %Out of bounds processing
                for j=1:llDim
                    if Ydec(j)<llDimMin(j)
                        Ydec(j)=llDimMin(j);

                    elseif Ydec(j)>llDimMax(j)
                        Ydec(j)=llDimMax(j);
                    end
                end
                llPop=[llPop;Ydec];
            end
        end
        %Handling xl0 out of bounds
        for j=1:llDim
            if xl0(j)<llDimMin(j)
                xl0(j)=llDimMin(j);
            elseif xl0(j)>llDimMax(j)
                xl0(j)=llDimMax(j);
            end
        end

        llPop=[llPop;xl0];

    end

    %llevaluate
    [ll_obj, ~, ll_con]=llevaluate(llPop, Problem, xu);
    llFunctionEvaluations = llFunctionEvaluations+llPopSize;

    %L_fitness
    if isempty(ll_con)
        ll_fitness=ll_obj;

    else
        PopCon   = sum(max(0,ll_con),2);
        Feasible = PopCon <= 0;
        ll_fitness  = Feasible.*ll_obj + ~Feasible.*(PopCon+1e10);
    end

    %Record the optimal value trajectory of the population
    [~,index]=min(ll_fitness);
    current_best_f=ll_obj(index,:);
    history_best_f=current_best_f;

    v0 = var(llPop);
    alpha = 1;
    GEN=0;

    %%Iterate
    while alpha>llStoppingCriteria && GEN<=llMaxGens

        %TournamentSelection
        MatingPool = TournamentSelection(2,llPopSize,ll_fitness);
        %produce offspring
        Offspring  = BLSAEAGA(llPop(MatingPool',:),{proC,disC,proM,disM},llDimMin, llDimMax);

        %llevaluate
        [off_ll_obj, ~, off_ll_con]=llevaluate(Offspring, Problem, xu);
        llFunctionEvaluations = llFunctionEvaluations+llPopSize;
        %L_fitness
        if isempty(off_ll_con)
            off_ll_fitness=off_ll_obj;

        else
            PopCon   = sum(max(0,off_ll_con),2);
            Feasible = PopCon <= 0;
            off_ll_fitness  = Feasible.*off_ll_obj + ~Feasible.*(PopCon+1e10);
        end

        %Merge parent and child populations
        llPop = [llPop;Offspring];ll_fitness=[ll_fitness;off_ll_fitness];  ll_obj=[ll_obj;off_ll_obj];
        [~,rank]   = sort(ll_fitness);
        llPop = llPop(rank(1:llPopSize),:);  ll_obj = ll_obj(rank(1:llPopSize),:);  ll_fitness = ll_fitness(rank(1:llPopSize),:);

        current_best_xl=llPop(1,:);

        %SQP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ub = Problem.upper(Problem.DU+1:end);
        lb = Problem.lower(Problem.DU+1:end);

        options = optimoptions('fmincon','Algorithm','sqp','Display','none');
        % Construct the Lower quadratic approximations of objective function and linear approximations of constraints
        approx.function       = QuadApprox(ll_obj, llPop);
        approx.equalityConstr = [];
        if size(off_ll_con,2) ~= 0
            for i = 1 : size(off_ll_con,2)
                approx.inequalityConstr{i} = QuadApprox(off_ll_con(:,i), llPop);
            end
        else
            approx.inequalityConstr = [];
        end
        [x,~,~,output] = fmincon(@(x) approximatedFunction(x,approx.function),current_best_xl,[],[],[],[],lb,ub,@(x) approximatedConstraints(x,approx.equalityConstr,approx.inequalityConstr),options);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        L_current_best_xl=x;
        llFunctionEvaluations1=output.funcCount;

        llFunctionEvaluations=llFunctionEvaluations+llFunctionEvaluations1;

        [L_off_ll_obj, ~, L_off_ll_con]=llevaluate(L_current_best_xl, Problem, xu);
        if isempty(L_off_ll_con)
            L_off_ll_fitness=L_off_ll_obj;

        else
            PopCon   = sum(max(0,L_off_ll_con),2);
            Feasible = PopCon <= 0;
            L_off_ll_fitness  = Feasible.*L_off_ll_obj + ~Feasible.*(PopCon+1e10);
        end

        if ll_fitness(llPopSize,:)>L_off_ll_fitness
            llPop(llPopSize,:)=L_current_best_xl;
            ll_obj(llPopSize,:)=L_off_ll_obj;
            ll_fitness(llPopSize,:)=L_off_ll_fitness;
        end

        alpha = sum(var(llPop)./v0);
        if alpha>1
            alpha = 1;
        end

        GEN=GEN+1;

        [~,index]=min(ll_fitness);
        current_best_f=ll_obj(index,:);
        history_best_f=[history_best_f current_best_f];


        if GEN>9
            if  abs(history_best_f(end)-history_best_f(end-5))<1e-5
                break;
            end

        end

    end

    val=min(ll_fitness);
    value=abs(ll_fitness-val);
    k=find(value<1e-9);
    llPop=llPop(k,:);    ll_obj=ll_obj(k,:);
    [n,~]=size(llPop);
    A=[];
    for i=1:n
        [ul_obj, ~, ul_con]=ulevaluate(xu, llPop(i,:), Problem);
        A=[A ul_obj];
    end

    [~,index]=min(A);
    best_xl=llPop(index,:);
    best_f=ll_obj(index,:);

    [~, ~, best_con]=llevaluate(best_xl, Problem, xu);

    if isempty(best_con)
        best_con=0;
    end
end

function approxFunctionValue = approximatedFunction(pop, parameters)
    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop';
end

%approximatedConstraints
function [c, ceq] = approximatedConstraints(pop, parametersEqualityConstr, parametersInequalityConstr)
    if ~isempty(parametersEqualityConstr)%等式约束
        for i = 1 : length(parametersEqualityConstr)
            ceq(i) = parametersEqualityConstr{i}.constant + pop*parametersEqualityConstr{i}.linear + pop*parametersEqualityConstr{i}.sqmatrix*pop';
        end
    else
        ceq = [];
    end

    if ~isempty(parametersInequalityConstr)%不等式约束
        for i = 1 : length(parametersInequalityConstr)
            c(i) = parametersInequalityConstr{i}.constant + pop*parametersInequalityConstr{i}.linear + pop*parametersInequalityConstr{i}.sqmatrix*pop';
        end
    else
        c = [];
    end
end

%OperatorGA - Crossover and mutation operators of genetic algorithm.
function Offspring = BLSAEAGA(Parent,Parameter,llDimMin, llDimMax)
    % Parameter setting
    if nargin > 1
        [proC,disC,proM,disM] = deal(Parameter{:});
    else
        [proC,disC,proM,disM] = deal(1,20,1,20);
    end

    Parent1 = Parent(1:floor(end/2),:);% get the two step reason
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);
  
   
    % Genetic operators for real encoding  % Simulated binary crossover
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                         (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    % Polynomial mutation
    Lower = repmat(llDimMin,2*N,1);
    Upper = repmat(llDimMax,2*N,1);
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end