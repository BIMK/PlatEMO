function Archive = UpdateArchive(N,combinePopulation)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen

    % Remove the dominated solutions
    Archive = combinePopulation(NDSort(combinePopulation.objs,1)==1);
    
    % Update the Archive outPopulation, if the number of solutions is larger than the population size
    if length(Archive) > N
        PopObj = Archive.objs;
        
        Choose = false(1,size(PopObj,1));       
        % Select the extreme solutions
        [~,extreme]     = min(pdist2(PopObj,eye(size(PopObj,2)),'cosine'),[],1);
        Choose(extreme) = true;
        
        %% Lp-norm-distances between each two solutions
        LpNormD = pdist2(PopObj,PopObj,'minkowski',0.5);
        while sum(Choose) < N
            Remain   = find(~Choose);
            [~, rho] = max(min(LpNormD(Remain,Choose),[],2));
            Choose(Remain(rho)) = true;
        end
        Archive = Archive(Choose);
    end
end