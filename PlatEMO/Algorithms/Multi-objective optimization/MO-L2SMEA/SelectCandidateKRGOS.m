function newOffspring = SelectCandidateKRGOS(subproblemList,nLinear,D,candidates)
% A surrogate-assisted offspring generation method for expensive multi-objective optimization problems

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Get surrogate fitness
    PopObj  = candidates(:,end-nLinear+1:2*nLinear);
    
    %% Non-dominated sorting
    FNO = NDSort(PopObj,size(PopObj,1));
    NDOffspring = candidates(FNO==1,:);
    
    %% Candidates selection process of each surrogate model
    newOffspring = [];
    for k = 1 : nLinear
        % Delete candidates closing to existing solutions
        theta   = 0.07 *(subproblemList{k}.ub-subproblemList{k}.lb)*nLinear;
        dist    = pdist2(NDOffspring(:,k),subproblemList{k}.trainDec);
        minDist = min(dist,[],2);
        remain  = minDist >= theta;
        remainDec = NDOffspring(remain,k);
        remainFit = NDOffspring(remain,nLinear+k);
        remainStd = -NDOffspring(remain,2*nLinear+k);
        % Delete solutions closing to new candidate solutions
        if sum(remain) > 2
            [frontNo,~] = NDSort([remainStd,remainFit],1);
            select      = 1:1:sum(remain);
            while any(frontNo>1) && length(select) > 2
                delIdx  = frontNo>1;
                select(delIdx) = [];
                [frontNo,~] = NDSort([remainStd(select),remainFit(select)],1);
            end
            finalDec = remainDec(select);
            finalObj = remainFit(select);
            n = length(select);
            popDec = repmat(subproblemList{k}.start,n,1) + repmat(finalDec,1,D).*repmat(subproblemList{k}.direct,n,1);
            newOffspring = [newOffspring;popDec,finalObj];
        elseif sum(remain) >= 1
            finalDec = remainDec;
            finalObj = remainFit;
            n = sum(remain);
            popDec = repmat(subproblemList{k}.start,n,1) + repmat(finalDec,1,D).*repmat(subproblemList{k}.direct,n,1);
            newOffspring = [newOffspring;popDec,finalObj];
        end
    end
    if isempty(newOffspring)
        [MR,MI] = min(NDOffspring(:,nLinear+1:2*nLinear));
        [~,LI]  = min(MR);
        if size(NDOffspring,1) == 1
            finalDec = NDOffspring(1,MI);
            finalObj = NDOffspring(1,nLinear+MI);
            popDec = repmat(subproblemList{MI}.start,1,1) + repmat(finalDec,1,D).*repmat(subproblemList{MI}.direct,1,1);
        else
            finalDec = NDOffspring(MI(LI),LI);
            finalObj = NDOffspring(MI(LI),nLinear+LI);
            popDec   = repmat(subproblemList{LI}.start,1,1) + repmat(finalDec,1,D).*repmat(subproblemList{LI}.direct,1,1);
        end
        newOffspring = [popDec,finalObj];
    end
end