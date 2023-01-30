function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Problem,Population,N,type)
% Environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if type
        %% Optimization of HDRVs
        DRI = zeros(length(Population),1);
        for i = 1 : length(Population)
            PopX = Problem.Perturb(Population(i).dec);
            [FrontNo,MaxFNo] = NDSort(PopX.objs,inf);
            for j = 1 : MaxFNo
                DRI(i) = DRI(i) + numel(find(FrontNo==j))*(j-1);
            end
            DRI(i) = DRI(i)/length(PopX);
        end
        DRI    = (DRI-min(DRI))./(max(DRI)-min(DRI));
        PopObj = Population.objs + repmat(1./(1-DRI+1e-6),1,Problem.M);
        % Non-dominated sorting
        [FrontNo,MaxFNo] = NDSort(PopObj,N);
        Next = FrontNo < MaxFNo;
        % DRI based selection
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(DRI(Last));
        Next(Last(Rank(1:N-sum(Next)))) = true;
        % Population for next generation
        Population = Population(Next);
        FrontNo    = FrontNo(Next);
        CrowdDis   = DRI(Next);
    else
        %% Optimization of LDRVs
        % Non-dominated sorting
        [FrontNo,MaxFNo] = NDSort(Population.objs,N);
        Next = FrontNo < MaxFNo;
        % Crowding distance based selection
        CrowdDis = CrowdingDistance(Population.objs,FrontNo);
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:N-sum(Next)))) = true;
        % Population for next generation
        Population = Population(Next);
        FrontNo    = FrontNo(Next);
        CrowdDis   = CrowdDis(Next);
    end
end