function [Population, Samples, winning_weights] = Operator(Problem, Population, W, winning_weights)    

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias (email: lrcf@cin.ufpe.br)

    % Generate all candidate offsprings
    MatingPool   = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2))';
    CandidateDec = Population.decs;
    OffspringDec = OperatorDE(Problem,CandidateDec(MatingPool,:), CandidateDec(randi(Problem.N,1,Problem.N),:) , CandidateDec(randi(Problem.N,1,Problem.N),:));    

	% Normalize offspring for SOM mapping
	Normalized_OffspringDec = rescale(OffspringDec,'InputMin',Problem.lower,'InputMax',Problem.upper);
	
	% Distance between each solution to a neuron
    Distance = pdist2(Normalized_OffspringDec,W(:,1:Problem.D)); 
    [~,rank] = sort(Distance,2);    
    
	% Estimate the offspring objective values using SOM
	Offspring_Labels = W(rank(:,1),Problem.D+1:Problem.D+Problem.M);

	%% Assign one to winning nodes
    winning_weights(rank(:,1)) = true;	
	
	%% Selection of survivors by ranking
	objs = [Population.objs;Offspring_Labels];
	cons = [Population.cons;zeros(size(Offspring_Labels,1),size(Population.cons,2))];
    
	%% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(objs,cons,Problem.N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:Problem.N-sum(Next)))) = true;
    
    %% Population for next generation
    out = Next(1:Problem.N);
    in  = Next(Problem.N+1:end);

	%% Selection
    if sum(in) >= 1
        Offspring = Problem.Evaluation(OffspringDec(in,:));
        Samples   = Offspring;
        Population(~out) = Offspring;
    else
        Samples = Population;
    end
end