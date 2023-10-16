function eliteIndiv = llSearch(Problem,ulPopDec,llPopDec)
% Obtain the upper member corresponds to the best lower member

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Lower level population initialization
    llPopDec     = [llPopDec;unifrnd(repmat(Problem.lower(Problem.DU+1:end),Problem.N-size(llPopDec,1),1),repmat(Problem.upper(Problem.DU+1:end),Problem.N-size(llPopDec,1),1))];
    llPopulation = Problem.EvaluationLower([repmat(ulPopDec,Problem.N,1),llPopDec]);
    FElower      = 0;
    
    %% Optimization
    while FElower < Problem.maxFElower
        % Select parents and generate offspring
        MatingPool  = TournamentSelection(2,3,CalFitness(Problem.C,llPopulation));
        ParentDec   = llPopulation(MatingPool).decs;
        llOffDec    = OperatorPCX(ParentDec(:,Problem.DU+1:end),Problem.lower(Problem.DU+1:end),Problem.upper(Problem.DU+1:end));
        llOffspring = Problem.EvaluationLower([repmat(ulPopDec,size(llOffDec,1),1),llOffDec]);
        FElower     = FElower + length(llOffspring);
        % Select r members with better adaptability
        llPopulation = EnvironmentalSelection(Problem,llPopulation,llOffspring);
    end
    [~,best]   = min(CalFitness(Problem.C,llPopulation));
    eliteIndiv = llPopulation(best).dec(Problem.DU+1:end);
end