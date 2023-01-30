function Population = FinalSelection(Archive,W,ArcW,ArcSP)
% Select final robust solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the convergence performance of solutions
    genx = cellfun(@length,ArcW);
    A_SP = cellfun(@(s)sum(s-min(s)),ArcSP);
    SM   = mean(A_SP);
    GM   = min(30,mean(genx));
    S    = cellfun(@mean,ArcSP) + SM/GM*max(0,GM-genx);
    
    %% Select a solution for each weight vector
    Selected = zeros(1,size(W,1));
    while true
        % Find a solution for each weight vector
        ClosestW  = cellfun(@mode,ArcW);
        AssignedW = [];
        for i = find(Selected==0)
            current = find(ClosestW==i);
            if ~isempty(current)
                [~,best]    = min(S(current));
                Selected(i) = current(best);
                AssignedW   = [AssignedW,i];
            end
        end
        % Delete the assigned weight vectors
        ArcW = cellfun(@(s)s(~ismember(s,AssignedW)),ArcW,'UniformOutput',false);
        if isempty(AssignedW)
            break;
        end
    end
    Population = Archive(Selected(Selected~=0));
end