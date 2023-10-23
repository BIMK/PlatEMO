function [TArchive,Archive] = UpdateArchive(TArchive,Problem,popNew,Archive)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    newN     = size(popNew,1);
    D        = Problem.D;
    BU       = Problem.upper;
    BD       = Problem.lower;
    RemainFE = Problem.maxFE - Problem.FE;
    
    %% Update Archive
    if newN == 0
        return;
    elseif newN <= RemainFE
        % Evaluate all candidate solutions with true function
        Dec         = popNew(:,1:D).*repmat(BU-BD,newN,1) + repmat(BD,newN,1);
        NewSolution = Problem.Evaluation(Dec);
        TArchive    = [TArchive;popNew(:,1:D),NewSolution.objs];
        Archive     = [Archive, NewSolution];
    else
        % Evaluate top 'remain' candidate solutions with true function
        [~,idx]     = sort(popNew(:,D+1),'ascend');
        topRemain   = idx(1:RemainFE);
        Dec         = popNew(topRemain,1:D).*repmat(BU-BD,RemainFE,1) + repmat(BD,RemainFE,1);
        NewSolution = Problem.Evaluation(Dec);
        TArchive    = [TArchive;popNew(topRemain,1:D),NewSolution.ojs];
        Archive     = [Archive, NewSolution];
    end
end