function Arc = WeightOptimization(Problem,G2,Population,wD,N)
% The second step of LSMOF, which aims to search the PF according to the
% bi-direction weight vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

	%% Choose NF solutions as the reference solutions
	Reference    = max(Population.objs,[],1);
    [RefPop,~,~] = EnvironmentalSelection(Population,wD);
    
    %% Calculate the reference directions
	Direction = [sum((RefPop.decs-repmat(Problem.lower,wD,1)).^2,2).^(0.5);sum((repmat(Problem.upper,wD,1)-RefPop.decs).^2,2).^(0.5)];
	Direct    = [(RefPop.decs-repmat(Problem.lower,wD,1));(repmat(Problem.upper,wD,1)-RefPop.decs)]./repmat(Direction,1,Problem.D);
	wmax      = sum((Problem.upper-Problem.lower).^2)^(0.5)*0.5;
    
    %% Optimize the weight variables by DE
	w0 = rand(N,2*wD).*wmax;                                    % Initialize the population
    [fitness,PopNew] = fitfunc(Problem,w0,Direct,Reference);	% Calculate the fitness and store the solutions
	Arc = PopNew(NDSort(PopNew.objs,1)==1);
	pCR = 0.2;
    beta_min=0.2;   % Lower Bound of Scaling Factor
    beta_max=0.8;   % Upper Bound of Scaling Factor
    empty_individual.Position=[];
    empty_individual.Cost=[];
    pop=repmat(empty_individual,N,1);
    for i = 1 : N
        pop(i).Position = w0(i,:);
        pop(i).Cost = fitness(i);
    end
    temp = [];
    for it = 1 : G2
        for i = 1 : N
            x = pop(i).Position;
            A = randperm(N);
            A(A==i) = [];
            a = A(1); b = A(2); c = A(3);
            % Mutation  %beta=unifrnd(beta_min,beta_max);
            beta = unifrnd(beta_min,beta_max,[1 2*wD]);
            y = pop(a).Position + beta.*(pop(b).Position - pop(c).Position);
            y = min(max(y,0),wmax);
            % Crossover
            z = zeros(size(x));
            j0=randi([1 numel(x)]);
            for j=1:numel(x)
                if j==j0 || rand<=pCR
                    z(j) = y(j);
                else
                    z(j) = x(j);
                end
            end
            NewSol.Position = z;
            [fit,PopNew] = fitfunc(Problem,z,Direct,Reference);
            temp = [temp,PopNew];
            temp = temp(NDSort(temp.objs,1)==1);
            NewSol.Cost = fit;
            if NewSol.Cost < pop(i).Cost
                pop(i)=NewSol;
            end
        end
    end
    %Update and store the non-dominated solutions
    Arc = [Arc,temp];
    if length(Arc) > Problem.N
        [frontNo,~] = NDSort(Arc.objs,1);
        Arc = Arc(frontNo==1);
    end
end

function [Obj,OffSpring] = fitfunc(Problem,w0,direct,Reference)
    [SubN,WD] = size(w0); 
    WD        = WD/2;
    Obj   	  = zeros(SubN,1);
    OffSpring = [];
    for i = 1 : SubN 
        PopDec    = [repmat(w0(i,1:WD)',1,Problem.D).*direct(1:WD,:)+repmat(Problem.lower,WD,1);
                     repmat(Problem.upper,WD,1) - repmat(w0(i,WD+1:end)',1,Problem.D).*direct(WD+1:end,:)];
        OffWPop   = Problem.Evaluation(PopDec);
        OffSpring = [OffSpring,OffWPop];
        Obj(i)    = -HV(OffWPop,Reference);
    end
end