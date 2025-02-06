function [Population,Dec,Mask] = BinaryGroupOptimization(Problem,Population,Dec,Mask,r_eva,nGroup,REAL)
% Binary group optimization framework

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Choose nondominated solutions as the reference solutions
    if REAL
        [~,~,RefDec,~,~] = EnvironmentalSelection(Population,Dec,Mask,Problem.M+1);
    else
        
        RefDec = ones(1,Problem.D);
    end
    
    % Divide variables with grouping technique and optimize the binary vectors
    [Subcomponent,TPop,TDec,TMask] = Group(Problem,REAL,nGroup);
    eval = r_eva * Problem.maxFE - Problem.FE;
    for i = 1 : size(RefDec,1)
        
        % Initilization
        W      = randi([0 1],Problem.N,length(Subcomponent));
        W(1,:) = [1,zeros(1,length(Subcomponent)-1)];
        Pop    = fitfunc(Problem,W,Subcomponent,RefDec(i,:));
        [~,~,~,W_FrontNo,W_CrowdDis] = EnvironmentalSelection(Pop,repmat(RefDec(i,:),[Problem.N,1]),W,Problem.N);
        
        % Main loop
        for j = 1 : floor(eval/Problem.N/size(RefDec,1))
            MatingPool = TournamentSelection(2,Problem.N,W_FrontNo,-W_CrowdDis);
            [~,OffW]   = Operator(Problem,repmat(RefDec(i,:),[Problem.N,1]),W(MatingPool,:),REAL);
            OffPop     = fitfunc(Problem,OffW,Subcomponent,RefDec(i,:));
            [Pop,~,W,W_FrontNo,W_CrowdDis] = EnvironmentalSelection([Pop,OffPop],repmat(RefDec(i,:),[2*Problem.N,1]),[W;OffW],Problem.N);
        end
        
        % Update the population
        OffMask = zeros(size(W,1),Problem.D);
        OffDec  = OffMask;
        for j = 1 : size(W,1)
            Selected = Subcomponent(W(j,:)==1);
            for k = 1 : length(Selected)
                OffMask(j,Selected{k}) = 1;
            end
            if REAL
                OffDec(j,:) = Match(Dec,OffMask(j,:),Problem);
            else
                OffDec(j,:) = ones(1,Problem.D);
            end
        end
        Offspring  = Problem.Evaluation(OffDec.*OffMask);
        [Population,Dec,Mask,~,~] = EnvironmentalSelection([Population,Offspring,TPop],[Dec;OffDec;TDec],[Mask;OffMask;TMask],Problem.N);
        drawnow();
    end
end

function [Subcomponent,Pop,Dec,Mask] = Group(Problem,REAL,nGroup)
% Grouping with fitness
    Fitness = zeros(1,Problem.D);
    for i = 1 : 1+4*REAL
        if REAL
            Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
        else
            Dec = ones(Problem.D,Problem.D);
        end
        Mask = eye(Problem.D);
        Pop  = Problem.Evaluation(Dec.*Mask);
        Fitness = Fitness + NDSort([Pop.objs,Pop.cons],inf);
    end
    
    Subcomponent = cell(1,nGroup);
    gamma = ceil(Problem.D/nGroup);
    j     = 1;
    while sum(Fitness~=Inf) > gamma
        [~,I] = sort(Fitness);
        Subcomponent{j} = I(1:gamma);
        Fitness(Subcomponent{j}) = Inf;
        j = j + 1;
    end
    Subcomponent{j} = find(Fitness~=Inf);
end

function Offspring = fitfunc(Problem,W,Subcomponent,RefPop)
% Calculate the objective values based on the weights and reference population

    [SubN,~]  = size(W);
    Offspring = [];
    for i = 1 : SubN
        Mask     = zeros(1,Problem.D);
        Selected = Subcomponent(W(i,:)==1);
        for j = 1 : length(Selected)
            Mask(Selected{j}) = 1;
        end
        
        PopDec    = RefPop.*Mask;
        OffWPop   = Problem.Evaluation(PopDec);
        Offspring = [Offspring,OffWPop];
    end
end