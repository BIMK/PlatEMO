classdef BLSAEA < ALGORITHM
% <2024> <multi> <real> <constrained/none> <bilevel>
% Bi-level surrogate modelling based evolutionary algorithm

%------------------------------- Reference --------------------------------
% H. Jiang, K. Qiu, Y. Tian, X. Zhang, and Y. Jin. Efficient surrogate
% modeling method for evolutionary algorithm to solve bilevel optimization
% problems. IEEE Transactions on Cybernetics, 2024, 54(7): 4335-4347
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)

            ulPopSize=floor(Problem.N/2);                   % Size of UL population
            ulDim=Problem.DU;                               % Number of UL dimensions
            llPopSize=Problem.N-ulPopSize;                  % Size of LL population
            llMaxGens=2000;                                 % Maximum number of generations allowed at LL
            llDim=Problem.DL;                               % Number of LL dimensions

            ulDimMin=Problem.lower(1:Problem.DU);           % Minimum bound accross UL dimensions
            ulDimMax=Problem.upper(1:Problem.DU);           % Maximum bound accross UL dimensions

            llDimMin=Problem.lower(Problem.DU+1:end);    	% Minimum bound accross LL dimensions
            llDimMax=Problem.upper(Problem.DU+1:end);    	% Maximum bound accross LL dimensions

            llStoppingCriteria = 1e-5;

            %% Ulsearch
            proC=1;disC=20;proM=1;disM=20;

            ulFunctionEvaluation = 0;
            llFunctionEvaluation = 0;

            %Upper population initialization (LHS sampling)
            Pu = UniformPoint(ulPopSize,ulDim,'Latin');
            ulPop = repmat(ulDimMax-ulDimMin,ulPopSize,1).*Pu+repmat(ulDimMin,ulPopSize,1);

            %For each xu, find the optimal lower-level solution
            for i=1:ulPopSize
                [llPop(i,:), llFunctionValue(i,:),l_con(i,:),llFunctionEvaluation,ulFunctionEvaluation]=llSearch([],0,[],[],Problem,llDim, llDimMin, llDimMax, ulPop(i,:),llPopSize,llStoppingCriteria,proC,disC,proM,disM,llFunctionEvaluation,llMaxGens,ulFunctionEvaluation);
            end

            %ulevaluate
            [ul_obj, ~, ul_con]=ulevaluate(ulPop, llPop, Problem);
            ulFunctionEvaluation = ulFunctionEvaluation+ulPopSize;
            if isempty(ul_con)
                [n,~]=size(ul_obj);
                ul_con=zeros(n,1);
            end

            %Fitness
            if isempty(ul_con) && all(l_con==0)
                ul_fitness=ul_obj;
            else
                PopCon   = sum(max(0,ul_con),2)+sum(max(0,l_con),2);
                Feasible = PopCon <= 0;
                ul_fitness  = Feasible.*ul_obj + ~Feasible.*(PopCon+1e10);
            end

            %Archive
            Archive_xu=ulPop; Archive_xl=llPop; Archive_F=ul_obj; Archive_f=llFunctionValue;

            GEN=0;

            Population  = Problem.Evaluation([ulPop,llPop]);
            %Iterate
            while Algorithm.NotTerminated(Population)
                %Selection
                MatingPool1 = TournamentSelection(2,ulPopSize,ul_fitness);
                MatingPool2 = TournamentSelection(2,ulPopSize,ul_fitness);

                %Produce offspring
                Offspring1  = BLSAEAGA(ulPop(MatingPool1',:),{proC,disC,proM,disM},ulDimMin, ulDimMax);
                Offspring2  = BLSAEAGA(ulPop(MatingPool2',:),{proC,disC,proM,disM},ulDimMin, ulDimMax);

                Offspring  = [Offspring1 ; Offspring2];

                X=Archive_xu;
                Y=Archive_F;

                % srgtsKRGSetOptions
                options = srgtsKRGSetOptions(X,Y);

                % srgtsKRGFit
                switch func2str(options.FIT_Fn)
                    case 'dace_fit'
                        state.FIT_Fn    = options.FIT_Fn;
                        state.FIT_FnVal = NaN;
                        if isempty(options.KRG_LowerBound) % no optimization for theta
                            [surrogate.KRG_DACEModel, state.KRG_DACEPerf, state.FIT_FnVal] = dace_fit(...
                                options.P, ...
                                options.T, ...
                                options.KRG_RegressionModel, ...
                                options.KRG_CorrelationModel, ...
                                options.KRG_Theta0);
                        else
                            [surrogate.KRG_DACEModel, state.KRG_DACEPerf, state.FIT_FnVal] = dace_fit(...
                                options.P, ...
                                options.T, ...
                                options.KRG_RegressionModel, ...
                                options.KRG_CorrelationModel, ...
                                options.KRG_Theta0, ...
                                options.KRG_LowerBound, ...
                                options.KRG_UpperBound);
                        end

                    case 'srgtsXVFit'
                        [surrogate state] = srgtsXVFit(options);
                end

                [pre_off_ul_obj,PredVar] = dace_predictor(Offspring, surrogate.KRG_DACEModel);

                objs=[pre_off_ul_obj PredVar];
                [FrontNo,~] = NDSort(objs,ulPopSize*2);

                %Select according to non-dominated level, first frontier
                Next = FrontNo ==1;
                Num = find(Next==1);
                non_n=length(Num);

                %For non-dominated descendants, find the optimal lower-level solution
                [Archive_size,~]=size(Archive_F);

                for j=1:non_n

                    %Calculate the cosine of the included angle and find the closest point K
                    for e=1:Archive_size
                        c_o_s(e,1)=1-pdist([Offspring(Num(j),:);Archive_xu(e,:)],'cosine');
                    end
                    c_o_s=c_o_s(1:Archive_size);

                    %Establish the mapping relationship of xu-xl(i), local model
                    models = struct('surrogate', {}, 'state', {});
                    xl0=[];
                    for i = 1:llDim
                        % srgtsKRGSetOptions
                        options = srgtsKRGSetOptions(X,Y);

                        % srgtsKRGFit
                        switch func2str(options.FIT_Fn)
                            case 'dace_fit'
                                models(i).state.FIT_Fn    = options.FIT_Fn;
                                models(i).state.FIT_FnVal = NaN;
                                if isempty(options.KRG_LowerBound) % no optimization for theta
                                    [models(i).surrogate.KRG_DACEModel, models(i).state.KRG_DACEPerf, models(i).state.FIT_FnVal] = dace_fit(...
                                        options.P, ...
                                        options.T, ...
                                        options.KRG_RegressionModel, ...
                                        options.KRG_CorrelationModel, ...
                                        options.KRG_Theta0);
                                else
                                    [models(i).surrogate.KRG_DACEModel, models(i).state.KRG_DACEPerf, models(i).state.FIT_FnVal] = dace_fit(...
                                        options.P, ...
                                        options.T, ...
                                        options.KRG_RegressionModel, ...
                                        options.KRG_CorrelationModel, ...
                                        options.KRG_Theta0, ...
                                        options.KRG_LowerBound, ...
                                        options.KRG_UpperBound);
                                end

                            case 'srgtsXVFit'
                                [models(i).surrogate models(i).state] = srgtsXVFit(options);
                        end
                        % srgtsKRGPredictor
                        [xl(i),~] = dace_predictor(Offspring(Num(j),:), models(i).surrogate.KRG_DACEModel);
                        xl0=[xl0 xl(i)];

                    end
                    %As the initial solution for lower-level optimization
                    [ND_off_llPop(j,:), ND_off_llFunctionValue(j,:),ND_l_con(j,:),llFunctionEvaluation,ulFunctionEvaluation]=llSearch(xl0,1,ulPop,llPop,Problem,llDim, llDimMin, llDimMax, Offspring(Num(j),:),llPopSize,llStoppingCriteria,proC,disC,proM,disM,llFunctionEvaluation,llMaxGens,ulFunctionEvaluation);

                end

                %ulevaluate
                [ND_off_ul_obj, ~, ND_off_ul_con]=ulevaluate(Offspring(Num',:), ND_off_llPop, Problem);
                ulFunctionEvaluation = ulFunctionEvaluation+non_n;
                if isempty(ND_off_ul_con)
                    [n,~]=size(ND_off_ul_obj);
                    ND_off_ul_con=zeros(n,1);
                end

                %Fitness
                if isempty(ND_off_ul_con) && all(ND_l_con(1:non_n,:)==0)
                    ND_off_ul_fitness=ND_off_ul_obj;
                else
                    PopCon   = sum(max(0,ND_off_ul_con),2)+sum(max(0,ND_l_con(1:non_n,:)),2);
                    Feasible = PopCon <= 0;
                    ND_off_ul_fitness  = Feasible.*ND_off_ul_obj + ~Feasible.*(PopCon+1e10);
                end

                %Newly generated descendant archive
                Archive_xu=[Archive_xu;Offspring(Num',:)]; Archive_xl=[Archive_xl;ND_off_llPop(1:non_n,:)]; Archive_F=[Archive_F;ND_off_ul_obj(1:non_n,:)]; Archive_f=[Archive_f;ND_off_llFunctionValue(1:non_n,:)];

                %Update archive
                A=round(Archive_xu, -9);
                [~,ia,~] = unique(A,'rows');

                Archive_xu = Archive_xu(ia,:);
                Archive_xl = Archive_xl(ia,:);
                Archive_F = Archive_F(ia,:);
                Archive_f = Archive_f(ia,:);
                [A_size,~]= size(Archive_F);
                if A_size>900
                    [~,F_rank]=sort(Archive_F);
                    Archive_xu=Archive_xu(F_rank(1:900),:);
                    Archive_xl=Archive_xl(F_rank(1:900),:);
                    Archive_F=Archive_F(F_rank(1:900),:);
                    Archive_f=Archive_f(F_rank(1:900),:);
                end

                %Merge parent and spring populations
                ulPop = [ulPop;Offspring(Num',:)];
                llPop = [llPop;ND_off_llPop(1:non_n,:)];
                ul_obj=[ul_obj;ND_off_ul_obj(1:non_n,:)];
                ul_fitness=[ul_fitness;ND_off_ul_fitness(1:non_n,:)];
                llFunctionValue=[llFunctionValue;
                    ND_off_llFunctionValue(1:non_n,:)];
                [~,rank]   = sort(ul_fitness);
                ulPop = ulPop(rank(1:ulPopSize),:);
                llPop = llPop(rank(1:ulPopSize),:);
                ul_obj = ul_obj(rank(1:ulPopSize),:);
                ul_fitness = ul_fitness(rank(1:ulPopSize),:);
                llFunctionValue = llFunctionValue(rank(1:ulPopSize),:);


                %Current optimal solution gradient
                current_best_xu=ulPop(1,:);

                [Archive_size,~]=size(Archive_F);

                %Calculate the cosine of the included angle and find the closest point K
                for e=1:Archive_size
                    c_o_s(e,1)=1-pdist([current_best_xu;Archive_xu(e,:)],'cosine');
                end
                c_o_s=c_o_s(1:Archive_size);

                %Establish the mapping relationship of xu-xl(i), local model
                models = struct('surrogate', {}, 'state', {});
                n_xl0=[];
                for i = 1:llDim
                    % srgtsKRGSetOptions
                    options = srgtsKRGSetOptions(X,Y);

                    % srgtsKRGFit
                    switch func2str(options.FIT_Fn)
                        case 'dace_fit'
                            models(i).state.FIT_Fn    = options.FIT_Fn;
                            models(i).state.FIT_FnVal = NaN;
                            if isempty(options.KRG_LowerBound) % no optimization for theta
                                [models(i).surrogate.KRG_DACEModel, models(i).state.KRG_DACEPerf, models(i).state.FIT_FnVal] = dace_fit(...
                                    options.P, ...
                                    options.T, ...
                                    options.KRG_RegressionModel, ...
                                    options.KRG_CorrelationModel, ...
                                    options.KRG_Theta0);
                            else
                                [models(i).surrogate.KRG_DACEModel, models(i).state.KRG_DACEPerf, models(i).state.FIT_FnVal] = dace_fit(...
                                    options.P, ...
                                    options.T, ...
                                    options.KRG_RegressionModel, ...
                                    options.KRG_CorrelationModel, ...
                                    options.KRG_Theta0, ...
                                    options.KRG_LowerBound, ...
                                    options.KRG_UpperBound);
                            end

                        case 'srgtsXVFit'
                            [models(i).surrogate models(i).state] = srgtsXVFit(options);
                    end
                    % srgtsKRGPredictor
                    [xl(i),~] = dace_predictor(Offspring(Num(j),:), models(i).surrogate.KRG_DACEModel);

                    n_xl0=[n_xl0 xl(i)];

                end

                %nested_LocalSearch
                lb = Problem.lower(1:Problem.DU);
                ub = Problem.upper(1:Problem.DU);

                options = optimoptions('fmincon','Algorithm','sqp','Display','none');
                Population = Problem.Evaluation([ulPop,llPop]);
                % Optimize the up-level objective using
                llPopulation = Problem.EvaluationLower([ulPop,llPop]);
                llPopCon = llPopulation.cons;
                llPopObj = llPopulation.objs;

                if sum(sum(isnan(llPopCon)))>0 || sum(sum(isnan(llPopObj)))>0
                    % llPopCon
                    nanColumnsCon = any(isnan(llPopCon), 1);
                    llPopCon(:, nanColumnsCon) = [];
                end
                l=size(llPopCon, 2);

                PopCon = Population.cons;
                PopObj = Population.objs;
                ulPopObj=PopObj(:,1);
                ulPopCon = PopCon(:,1:l);

                % Construct the Lower quadratic approximations of objective function and linear approximations of constraints
                approx.function       = QuadApprox(ulPopObj, ulPop);
                approx.equalityConstr = [];
                if size(ulPopCon,2) ~= 0
                    for i = 1 : size(ulPopCon,2)
                        approx.inequalityConstr{i} = QuadApprox(ulPopCon(:,i), ulPop);
                    end
                else
                    approx.inequalityConstr = [];
                end
                [x,~,~,output] = fmincon(@(x)nestedproblem(x, Problem,n_xl0,ulPop,llPop),current_best_xu,[],[],[],[],lb,ub,@(x) approximatedConstraints(x,approx.equalityConstr,approx.inequalityConstr),options);

                L_current_best_xu=x;
                ulFunctionEvaluation=ulFunctionEvaluation+output.funcCount;
                [L_xl, L_f,L_l_con,llFunctionEvaluation,ulFunctionEvaluation]=llSearch([],0,[],[],Problem,llDim, llDimMin, llDimMax, L_current_best_xu,llPopSize,llStoppingCriteria,proC,disC,proM,disM,llFunctionEvaluation,llMaxGens,ulFunctionEvaluation);
                [L_ul_obj, ~, L_ul_con,]=ulevaluate(L_current_best_xu, L_xl, Problem);
                % Fitness
                if isempty(L_ul_con) && all(L_l_con==0)
                    L_ul_fitness=L_ul_obj;
                else
                    PopCon   = sum(max(0,L_ul_con),2)+sum(max(0,L_l_con),2);
                    Feasible = PopCon <= 0;
                    L_ul_fitness  = Feasible.*L_ul_obj + ~Feasible.*(PopCon+1e10);
                end

                %Replacement population worst solution
                ulPop(ulPopSize,:)=L_current_best_xu;
                ul_obj(ulPopSize,:)=L_ul_obj;
                llPop(ulPopSize,:)=L_xl;
                llFunctionValue(ulPopSize,:)=L_f;
                ul_fitness(ulPopSize,:)=L_ul_fitness;

                GEN=GEN+1;

            end
        end
    end
end

function [ f,llFunctionEvaluations,xl] = nestedproblem(xa,Problem,n_xl0,ulPop,llPop)
    xu1 = xa(1:Problem.p);
    xu2 = xa(Problem.p+1:Problem.p+Problem.r);
    %SQP
    ub = Problem.upper(Problem.DU+1:end);
    lb = Problem.lower(Problem.DU+1:end); 

    options = optimoptions('fmincon','Algorithm','sqp','Display','none');

    llPopulation = Problem.EvaluationLower([ulPop,llPop]);
    llPopCon = llPopulation.cons;
    llPopObj = llPopulation.objs;
    if sum(sum(isnan(llPopCon)))>0 || sum(sum(isnan(llPopObj)))>0
        % Deal with llPopObj
        nanColumnsObj = any(isnan(llPopObj), 1);
        llPopObj(:, nanColumnsObj) = [];
        
        % Deal with llPopCon
        nanColumnsCon = any(isnan(llPopCon), 1);
        llPopCon(:, nanColumnsCon) = [];

    end
   % Construct the Lower quadratic approximations of objective function and linear approximations of constraints
    approx.function       = QuadApprox(llPopObj, llPop);
    approx.equalityConstr = [];
    if size(llPopCon,2) ~= 0
        for i = 1 : size(llPopCon,2)
            approx.inequalityConstr{i} = QuadApprox(llPopCon(:,i), llPop);
        end
    else
        approx.inequalityConstr = [];
    end
    [x,~,~,output] = fmincon(@(x) approximatedFunction(x,approx.function),n_xl0,[],[],[],[],lb,ub,@(x) approximatedConstraints(x,approx.equalityConstr,approx.inequalityConstr),options);
    xl=x;
    llFunctionEvaluations=output.funcCount;

    try
       xl1 = xl(1:Problem.q+Problem.s);
       xl2 = xl(Problem.q+Problem.s+1:Problem.q+Problem.s+Problem.r);
    catch
       xl1 = xl(1:Problem.q);
       xl2 = xl(Problem.q+1:Problem.q+Problem.r);
    end

    PopObj = Problem.CalObj([xu1 xu2 xl1 xl2]);
    f=PopObj(:,1);
end

function approxFunctionValue = approximatedFunction(pop, parameters)
    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop';
end

%approximatedConstraints
function [c, ceq] = approximatedConstraints(pop, parametersEqualityConstr, parametersInequalityConstr)
    if ~isempty(parametersEqualityConstr)
        for i = 1 : length(parametersEqualityConstr)
            ceq(i) = parametersEqualityConstr{i}.constant + pop*parametersEqualityConstr{i}.linear + pop*parametersEqualityConstr{i}.sqmatrix*pop';
        end
    else
        ceq = [];
    end

    if ~isempty(parametersInequalityConstr)
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