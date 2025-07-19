function [PCDec,PCObj,PCCon] = ArchiveUpdate(PopDec,PopObj,PopCon,N)
% Update archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select feasible solutions
    fIndex  = all(PopCon <= 0,2);
    FPopDec = PopDec(fIndex,:);
    FPopObj = PopObj(fIndex,:);
    FPopCon = PopCon(fIndex,:);
    if isempty(FPopDec)
        PCDec = [];
        PCObj = [];
        PCCon = [];
        return;
    else  
        NFPopDec = FPopDec(NDSort(FPopObj,1)==1,:);
        NFPopObj = FPopObj(NDSort(FPopObj,1)==1,:);
        NFPopCon = FPopCon(NDSort(FPopObj,1)==1,:);
        randind  = randperm(size(NFPopDec,1));
        PCDec    = NFPopDec(randind,:);
        PCObj    = NFPopObj(randind,:);
        PCObj2   = PCObj;
        PCCon    = NFPopCon(randind,:);
        nND      = size(NFPopDec,1);
        %% Population maintenance
        if size(PCDec,1) > N
            % Normalization
            fmax   = max(PCObj2,[],1);
            fmin   = min(PCObj2,[],1);
            PCObj2 = (PCObj2-repmat(fmin,nND,1))./repmat(fmax-fmin,nND,1);
            % Determine the radius of the niche
            d  = pdist2(PCObj2,PCObj2);
            d(logical(eye(length(d)))) = inf;
            sd = sort(d,2);
            r  = median(sd(:,min(size(PCObj2,2),size(sd,2))));
            R  = min(d./r,1);
            % Delete solution one by one
            while size(PCDec,1) > N
                [~,worst]      = max(1-prod(R,2));
                PCDec(worst,:) = [];
                PCObj(worst,:) = [];
                PCCon(worst,:) = [];
                R(worst,:)     = [];
                R(:,worst)     = [];
            end
        end
    end
end