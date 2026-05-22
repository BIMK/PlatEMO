classdef WOF_Utils
% WOF_Utils - Static class for all WOF utility and auxiliary functions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    methods(Static)
        %% 1. Create Global Dummy
        function GlobalDummy = create_global_dummy(gamma, xPrime, G, Global, populationSize, psi, optimiser)
            GlobalDummy       = {};
            GlobalDummy.lower = zeros(1,gamma);
            GlobalDummy.upper = ones(1,gamma).*2.0;
            if or(optimiser == 2, optimiser == 4)
                [uniW,GlobalDummy.N] = UniformPoint(populationSize,Global.M);
                GlobalDummy.uniW     = uniW;
            else
                GlobalDummy.N = populationSize;
            end
            GlobalDummy.xPrime      = xPrime;
            GlobalDummy.G           = G;
            GlobalDummy.psi         = psi;
            GlobalDummy.xPrimelower = Global.lower;
            GlobalDummy.xPrimeupper = Global.upper;
            GlobalDummy.isDummy     = true;
            GlobalDummy.Global      = Global;
        end
        
        %% 2. Eliminate Duplicates
        function Population = eliminate_duplicates(input)
            [~,ia,~]   = unique(input.objs,'rows');
            Population = input(ia);
        end
        
        %% 3. Fill Population
        function Population = fill_population(input, Problem)
            Population               = input;
            theCurrentPopulationSize = size(input,2);
            if theCurrentPopulationSize < Problem.N
                amountToFill = Problem.N-theCurrentPopulationSize;
                FrontNo      = NDSort(input.objs,inf);
                CrowdDis     = CrowdingDistance(input.objs,FrontNo);
                MatingPool   = TournamentSelection(2,amountToFill+1,FrontNo,-CrowdDis);
                Offspring    = OperatorGA(Problem,input(MatingPool));
                Population   = [Population,Offspring(1:amountToFill)];
            end
        end
        
        %% 4. Create Initial Weight Population
        function WeightPopulation = create_initial_weight_population(N, gamma, GlobalDummy)
            decs             = rand(N,gamma).*2.0;
            WeightPopulation = [];
            for i = 1:N
                solution         = WOF_WeightIndividual(decs(i,:),GlobalDummy);
                WeightPopulation = [WeightPopulation, solution];
            end
        end
        
        %% 5. Extract Population
        function W = extract_population(WeightPopulation, Problem, Population, G, psi, xPrime, q, methodToSelectxPrimeSolutions)
            % Step 1
            weightIndividualList = WOF_Utils.select_xprimes(WeightPopulation, q, methodToSelectxPrimeSolutions);
            calc                 = size(weightIndividualList,2)*size(Population,2);
            PopDec1              = ones(calc,Problem.D);
            count                = 1;
            for wi = 1 : size(weightIndividualList,2)
                weightIndividual = weightIndividualList(wi);
                weightVars       = weightIndividual.dec;
                for i = 1 : size(Population,2)
                    individualVars   = Population(i).dec;
                    x                = WOF_Utils.transformation_function_matrix_form(individualVars,weightVars(G),Problem.upper,Problem.lower, psi);
                    PopDec1(count,:) = x;
                    count            = count + 1;
                end
            end
            W1 = Problem.Evaluation(PopDec1);
            
            % Step 2
            PopDec2 = [];
            for wi = 1 : size(WeightPopulation,2)
                weightIndividual = WeightPopulation(wi);
                weightVars       = weightIndividual.dec;
                individualVars   = xPrime.dec;
                x                = 1:Problem.D;
                for j = 1 : Problem.D
                    x(j) = WOF_Utils.transformation_function(individualVars(j),weightVars(G(j)),Problem.upper(j),Problem.lower(j), psi);   
                end
                PopDec2 = [PopDec2;x]; 
            end
            W2 = Problem.Evaluation(PopDec2);
            W  = [W1,W2];
        end
        
        %% 6. Environmental Selection
        function [Population,FrontNo,CrowdDis] = environmental_selection(Population,N)
            [FrontNo,MaxFNo]     = NDSort(Population.objs,N);
            Next                 = false(1,length(FrontNo));
            Next(FrontNo<MaxFNo) = true;
            
            CrowdDis = CrowdingDistance(Population.objs,FrontNo);
            Last     = find(FrontNo==MaxFNo);
            [~,Rank] = sort(CrowdDis(Last),'descend');
            Popsize  = min(N,size(Population,2));
            Next(Last(Rank(1:Popsize-sum(Next)))) = true;
            
            Population = Population(Next);
            FrontNo    = FrontNo(Next);
            CrowdDis   = CrowdDis(Next);
        end
        
        %% 7. GA Half
        function Offspring = GA_half(Global,Parent)
            [proC,disC,proM,disM] = deal(1,20,1,20);
            Parent  = Parent.decs;
            Parent1 = Parent(1:floor(end/2),:);
            Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
            [N,D]   = size(Parent1);
            
            beta          = zeros(N,D);
            mu            = rand(N,D);
            beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
            beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
            beta          = beta.*(-1).^randi([0,1],N,D);
            beta(rand(N,D)<0.5)              = 1;
            beta(repmat(rand(N,1)>proC,1,D)) = 1;
            Offspring     = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;
            
            Lower = repmat(Global.lower,N,1);
            Upper = repmat(Global.upper,N,1);
            Site  = rand(N,D) < proM/D;
            mu    = rand(N,D);
            temp  = Site & mu<=0.5;
            
            Offspring       = min(max(Offspring,Lower),Upper);
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                              (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
            temp            = Site & mu>0.5; 
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                              (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
        end
        
        %% 8. Create Groups
        function [outIndexList,numberOfGroups] = create_groups(Problem,numberOfGroups,xPrime,method)
            switch method
                case 1 %linear grouping
                    varsPerGroup = floor(Problem.D/numberOfGroups);
                    outIndexList = [];
                    for i = 1 : numberOfGroups-1
                       outIndexList = [outIndexList, ones(1,varsPerGroup).*i];
                    end
                    outIndexList = [outIndexList, ones(1,Problem.D-size(outIndexList,2)).*numberOfGroups];
                case 2 %orderByValueGrouping
                    varsPerGroup = floor(Problem.D/numberOfGroups);
                    vars         = xPrime.dec;
                    [~,I]        = sort(vars);
                    outIndexList = ones(1,Problem.D);
                    for i = 1 : numberOfGroups-1
                       outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
                    end
                    outIndexList(I(((numberOfGroups-1)*varsPerGroup)+1:end)) = numberOfGroups;
                case 3 %random Grouping
                    varsPerGroup = floor(Problem.D/numberOfGroups);
                    outIndexList = [];
                    for i = 1 : numberOfGroups-1
                       outIndexList = [outIndexList, ones(1,varsPerGroup).*i];
                    end
                    outIndexList = [outIndexList, ones(1,Problem.D-size(outIndexList,2)).*numberOfGroups];
                    outIndexList = outIndexList(randperm(length(outIndexList)));
                case 4 %up or down groups
                    outIndexList = ones(1,Problem.D);
                    xPrimeVars   = xPrime.decs;
                    xPrimeObjs   = xPrime.objs;
                    for i = 1 : Problem.D
                        newSolVars    = xPrime.decs;
                        newSolVars(i) = xPrimeVars(i)*1.05;
                        newSol        = Problem.Evaluation(newSolVars);
                        newSolObjs    = newSol.objs;
                        if newSolObjs(1) < xPrimeObjs(1)
                            outIndexList(i) = 2;
                        end
                    end
                    numberOfGroups = 2;
            end
        end
        
        %% 9. Optimise By MOEAD
        function Population = optimise_by_MOEAD(GlobalDummy,Population,W,evaluations,isDummy)
            maximum = WOF_Utils.current_evaluations(GlobalDummy, isDummy) + evaluations;
            T       = max(ceil(GlobalDummy.N/10),2);
            B       = pdist2(W,W);
            [~,B]   = sort(B,2);
            B       = B(:,1:T);
            
            Z = min(Population.objs,[],1);
            g = zeros(GlobalDummy.N);
            for i = 1 : GlobalDummy.N
                g(i,:) = max(repmat(abs(Population(i).obj-Z),GlobalDummy.N,1)./W,[],2)';
            end
            [~,rank]  = sort(g,2);
            associate = zeros(1,GlobalDummy.N);
            for i = 1 : GlobalDummy.N
                x = find(~associate(rank(i,:)),1);
                associate(rank(i,x)) = i;
            end
            Population = Population(associate);
            
            while WOF_Utils.current_evaluations(GlobalDummy, isDummy) < maximum
                for i = 1 : GlobalDummy.N
                    P = B(i,randperm(size(B,2)));
                    if isDummy == true
                        NewDec    = WOF_Utils.GA_half(GlobalDummy, Population(P(1:2)));
                        Offspring = WOF_WeightIndividual(NewDec,GlobalDummy);
                    else
                        Offspring = OperatorGAhalf(GlobalDummy,Population(P(1:2)));
                    end
                    
                    Z    = min(Z,Offspring.obj);
                    type = 1;
                    switch type
                        case 1
                            normW   = sqrt(sum(W(P,:).^2,2));
                            normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
                            normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                            CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                            CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
                            g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                            g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                        case 2
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
                        case 3
                            Zmax  = max(Population.objs,[],1);
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                        case 4
                            g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                            g_new = max(repmat(abs(Offspring.obj-Z),T,1)./W(P,:),[],2);
                    end
                    Population(P(g_old>=g_new)) = Offspring;
                end
            end
        end
        
        %% 9.1 Helper for MOEAD
        function e = current_evaluations(GlobalDummy, isDummy)
            if isDummy == true  
                e = GlobalDummy.Global.FE;
            else
                e = GlobalDummy.FE;
            end
        end
        
        %% 10. Select xPrimes
        function weightIndList = select_xprimes(input,amount, method)
            inputSize = size(input,2);
            switch method 
                case 1 %largest Crowding Distance from first front
                    inFrontNo     = NDSort(input.objs,inf);
                    weightIndList = [];
                    i             = 1;
                    if inputSize < amount
                        weightIndList = input;
                    else
                        while size(weightIndList,2) < amount 
                            left          = amount - size(weightIndList,2);
                            theFront      = inFrontNo == i;
                            newPop        = input(theFront);
                            FrontNo       = NDSort(newPop.objs,inf);
                            CrowdDis      = CrowdingDistance(newPop.objs,FrontNo);
                            [~,I]         = sort(CrowdDis,'descend');
                            until         = min(left,size(newPop,2));
                            weightIndList = [weightIndList,newPop(I(1:until))];
                            i             = i + 1;
                        end
                    end
                case 2 %tournament selection by front and CD
                    FrontNo       = NDSort(input.objs,inf);
                    CrowdDis      = CrowdingDistance(input.objs,FrontNo);
                    weightIndList = input(TournamentSelection(2,amount,FrontNo,-CrowdDis));
                case 3 % first m+1 by reference lines + fill with random
                    objValues     = input.objs;
                    m             = size(objValues,2);
                    weightIndList = [];
                    for i = 1 : m
                        vec           = zeros(1,m);
                        vec(1,i)      = 1;
                        angles        = pdist2(vec,real(objValues),'cosine');
                        [~,minIndex]  = min(angles);
                        weightIndList = [weightIndList,input(minIndex)];
                    end
                    if size(weightIndList,2) < amount
                        vec           = ones(1,m);
                        angles        = pdist2(vec,real(objValues),'cosine');
                        [~,minIndex]  = min(angles);
                        weightIndList = [weightIndList,input(minIndex)];
                    end
                    while size(weightIndList,2) < amount
                        ind           = input(randi([1 inputSize],1,1));
                        weightIndList = [weightIndList,ind];
                    end
            end
        end
        
        %% 11. Transformation Function
        function value = transformation_function(xPrime,weight,maxVal,minVal,method)
            value = xPrime;
            switch method
                case 1 %multiplication
                    value = xPrime*weight;
                case 2 %p-value
                    pWert = 0.2;
                    value = xPrime+pWert*(weight-1.0)*(maxVal-minVal);
                case 3 %interval
                    if weight > 1.0
                        weight   = weight - 1.0;
                        interval = maxVal - xPrime;
                        value    = xPrime + weight * interval;
                    else
                        interval = xPrime - minVal;
                        value    = minVal + weight * interval;
                    end           
            end
            if value < minVal
               value = minVal;
            elseif value > maxVal
               value = maxVal;
            end
        end
        
        %% 12. Transformation Function Matrix Form
        function value = transformation_function_matrix_form(xPrime,weight,maxVal,minVal,method)
            value = xPrime;
            switch method
                case 1 %multiplication
                    value = xPrime*weight;
                case 2 %p-value
                    pWert = 0.2;
                    value = xPrime+pWert*(weight-1.0)*(maxVal-minVal);
                case 3 %interval
                    interval = xPrime - minVal;
                    value    = minVal + weight .* interval;
                    interval = maxVal - xPrime;
                    value(weight > 1.0) = xPrime(weight > 1.0) + (weight(weight > 1.0)-1.0) .* interval(weight > 1.0); 
            end
            if value < minVal
               value = minVal;
            elseif value > maxVal
               value = maxVal;
            end
        end
    end
end