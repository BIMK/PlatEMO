function [Population, PSI] = hpaEnvironmentalSelection(Population, Problem, V)
% The environmental selection of hpaEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen
    
    % Non-dominated sorting
    N = Problem.N;
    [FrontNo, MaxFNo] = NDSort(Population.objs, N);
    
    if numel(Population(FrontNo<=MaxFNo)) <= N
        SelI = FrontNo<=MaxFNo;
        Population = Population(SelI);
        
        Chooose = FrontNo==1;
        CC = Chooose + SelI;
        I = find(CC == 2);
        PSI = zeros(1, length(I));
        for j = 1: length(I)
            index = sum(SelI(1: I(j))); % the index of the solutions in the selected population
            PSI(j) = index;
        end
    else % If the number of selected solutions is larger than the population size       
        if MaxFNo <= 1
            
            [Population, PSI] = hpaNDSolutionsSelectionStrategy(Population(FrontNo==1), V, Problem);
                        
        else % retain some dominated solutions to improve diversity
            
            % Record the selected solutions
            SelectedSolutions = Population(FrontNo < MaxFNo);
            SelNum = numel(SelectedSolutions);
            % The solutions in the last front
            CanSolutions = Population(FrontNo==MaxFNo);
            CanNum = numel(CanSolutions);
            
            PopObj = [SelectedSolutions.objs; CanSolutions.objs];
            % Normalization
            Zmin   = min(PopObj, [], 1);
            Zmax   = max(PopObj, [], 1);
            PopObj = (PopObj-repmat(Zmin,size(PopObj,1),1))./repmat(Zmax-Zmin,size(PopObj,1),1);
            
            Population = [SelectedSolutions, CanSolutions];
            PSI = 1: SelNum; % The index of the porminent solutions
            Next = [true(1, SelNum), false(1, CanNum)];
            % Select remaining solutions based on angle
            angle = acos(1-pdist2(PopObj, PopObj, 'cosine'));
            while sum(Next) < N
                Select  = find(Next);
                Remain  = find(~Next);
                [~, rho] = max(min(angle(Remain, Select), [], 2));
                Next(Remain(rho)) = true;
            end
            Population = Population(Next);         
        end
    end
end