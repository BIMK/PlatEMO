function [Archive,RefPoint,Range, Ratio] = UpdateRefPoint(Archive,W,Range)
% Reference point adaption

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	%% Delete duplicated and dominated solutions
    Archive = unique(Archive(NDSort(Archive,1)==1,:),'rows');
    %unique(A,'rows') gets the matrix formed by the different vectors of matrix A
    NA      = size(Archive,1);
    NW      = size(W,1);
    
	%% Update the ideal point
    if ~isempty(Range)
        Range(1,:) = min([Range(1,:);Archive],[],1);
    elseif ~isempty(Archive)
        Range = [min(Archive,[],1);max(Archive,[],1)];%2*M
    end
    
    %% Update archive and reference points
    if size(Archive,1) <= 1
        RefPoint = W;
        Ratio=0;
    else
        %% Find contributing solutions and valid weight vectors
        tArchive = Archive - repmat(Range(1,:),NA,1);
        W        = W.*repmat(Range(2,:)-Range(1,:),NW,1);
        Distance      = CalDistance(tArchive,W); %The (1-cosine) value of the individual and the reference point W in tArchive, NA*NW
        [~,nearestP]  = min(Distance,[],1);
        ContributingS = unique(nearestP);
        [~,nearestW]  = min(Distance,[],2);
        ValidW        = unique(nearestW(ContributingS));%The reference point with the smallest Angle to the contributing solution

        %% Update archive
        Choose = ismember(1:NA,ContributingS);%1*NA Logical vector
        Cosine = 1 - pdist2(tArchive,tArchive,'cosine');%NA*NA
        Cosine(logical(eye(size(Cosine,1)))) = 0;% The Cosine matrix has a diagonal of 0
        while sum(Choose) < min(3*NW,size(tArchive,1))
            unSelected = find(~Choose);
            [~,x]      = min(max(Cosine(~Choose,Choose),[],2));
            % Calculate the minimum included Angle between all other solutions and each contributing solution, and add the remaining solutions with the maximum of the minimum included Angle into the contributing solution
            Choose(unSelected(x)) = true;
        end
        Archive  = Archive(Choose,:);
        tArchive = tArchive(Choose,:);

        %% Update reference points
        RefPoint = [W(ValidW,:);tArchive];
        Choose   = [true(1,length(ValidW)),false(1,size(tArchive,1))];
        Cosine   = 1 - pdist2(RefPoint,RefPoint,'cosine');
        Cosine(logical(eye(size(Cosine,1)))) = 0;
        while sum(Choose) < min(NW,size(RefPoint,1))
            Selected = find(~Choose);
            [~,x]    = min(max(Cosine(~Choose,Choose),[],2));
            Choose(Selected(x)) = true;
        end
        RefPoint = RefPoint(Choose,:);
        Ratio=length(ValidW)/NW;
    end 
end