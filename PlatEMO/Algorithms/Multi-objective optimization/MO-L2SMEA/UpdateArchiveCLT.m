function [Archive,currFEs,tArc,Population,MaxObj,MinObj] = UpdateArchiveCLT(Archive,Problem,popNew,remain,currFEs,tArc,Population,MaxObj,MinObj)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    newN = size(popNew,1);
    D    = Problem.D;
    BU   = Problem.upper;
    BD   = Problem.lower;
    
    %% Update Archive
    if newN == 0 || remain <= 0
        return;
    elseif newN <= remain
        % Evaluate all candidate solutions with true function
        Dec        = popNew(:,1:D).*repmat(BU-BD,newN,1) + repmat(BD,newN,1);
        NEWOFF     = Problem.Evaluation(Dec);
        Population = [Population,NEWOFF];
        DECS       = (NEWOFF.decs-repmat(Problem.lower,length(NEWOFF),1))./repmat(Problem.upper-Problem.lower,length(NEWOFF),1);
        MinObj     = min([MinObj;NEWOFF.objs],[],1);
        MaxObj     = max([MaxObj;NEWOFF.objs],[],1);
        tArc       = [tArc;DECS, NEWOFF.objs];
        Archive    = [Archive;DECS,NEWOFF.objs];
        currFEs    = currFEs + newN;
    else
        % Evaluate top 'remain' candidate solutions with true function
        [~,idx]    = sort(popNew(:,D+1),'ascend');
        topRemain  = idx(1:remain);
        Dec        = popNew(topRemain,1:D).*repmat(BU-BD,remain,1) + repmat(BD,remain,1);
        NEWOFF     = Problem.Evaluation(Dec);
        Population = [Population,NEWOFF];
        DECS       = (NEWOFF.decs-repmat(Problem.lower,length(NEWOFF),1))./repmat(Problem.upper-Problem.lower,length(NEWOFF),1);
        MinObj     = min([MinObj;NEWOFF.objs],[],1);
        MaxObj     = max([MaxObj;NEWOFF.objs],[],1);
        tArc       = [tArc;DECS, NEWOFF.objs];
        Archive    = [Archive;DECS,NEWOFF.objs];
        currFEs    = currFEs + remain;
    end
end