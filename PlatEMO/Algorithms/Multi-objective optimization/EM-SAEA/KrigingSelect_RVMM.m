function PopNew = KrigingSelect_RVMM(PopDec,PopObj,MSE,V,V0,NumV1,delta,mu_,theta,per,PopDec1,PopObj1,MSE1,V1,A2,A2Obj)
% Kriging selection in K-RVEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % nondominated sorting
    [c,ia,ic] = unique(PopObj1,'rows','stable');
    PopObj1   = PopObj1(ia,:);
    PopDec1   = PopDec1(ia,:);
    num       = 1;
    
    % flagMu denotes the number of solutions for real function evaluations
    flagMu = mu_;

    [c,ia,ic] = unique(PopObj,'rows','stable');
    PopObj    = PopObj(ia,:);
    PopDec    = PopDec(ia,:);
    Noid      = NDSort(PopObj,1)==1;
    PopObj    = PopObj(Noid,:);
    PopDec    = PopDec(Noid,:);
    
    Noid    = NDSort(PopObj1,1)==1;
    PopObj1 = PopObj1(Noid,:);
    PopDec1 = PopDec1(Noid,:);
    
    wholeobj = [A2Obj;PopObj1;PopObj];
    zmin     = min(wholeobj,[],1);
    Zmin     = min(A2Obj,[],1);
    zmin1    = min(wholeobj,[],1);

    scale = (max(A2Obj,[],1)-min(A2Obj,[],1));
    scale(:,scale==0) = 10^(-6);
    cd    = [];
    for i = 1 : size(PopObj1,1)
        mu    = PopObj1(i,:);
        R1{i} = mu;
        zmin1 = min([zmin1;R1{i}],[],1);
    end
    for i = 1 : size(PopObj1,1)
        d = pdist2((R1{i}-zmin1),(A2Obj - zmin1),'euclidean');
        A = R1{i};
        [cd_,id] = sort(d,2);
        id_ = id(:,1);
        k   = any(A<A2Obj(id_,:),2) - any(A>A2Obj(id_,:),2);
        cd_(k~=1,1) = 0;
        cd(i,:)     = sum(cd_(:,1));
    end
    if flagMu == 1
        [~,cbest] = max(cd,[],1);
        cid       = find(cd == cd(cbest,:));
        try
            cbest = cid(randperm(size(cid,1),1),:);
        catch e
            cbest = [];
        end
    elseif flagMu == 5
        [~,id] = sort(cd,1,'descend');
        cbest  = id(1:min(flagMu,size(cd,1)));
    end
    
    % diversity
    dbest = [];
    zmin  = min([zmin;zmin1],[],1);
    if ~isempty(PopObj)
        for i = 1 : size(PopObj,1)
            mu   = PopObj(i,:);
            R{i} = repmat(mu,num,1);
            zmin = min([zmin;R{i}],[],1);
        end
        allR = cell2mat(R(1:end)');
        if RVMM_IGD((max(allR-zmin,0))./scale,max( A2Obj- zmin,0)./scale) < RVMM_IGD((max(allR-zmin,0))./scale,max(A2Obj - Zmin,0)./scale)
            Angle = acos(1-pdist2((max(allR-zmin,0))./scale,max(A2Obj - zmin,0)./scale,'cosine'));
        else
            Angle = acos(1-pdist2((max(allR-zmin,0))./scale,max(A2Obj - Zmin,0)./scale,'cosine'));
        end
        [angle,~] = min(Angle,[],2);
        temp      = reshape(angle,num,length(R));
        dd        = mean(temp,1);
    
        dd = dd';
        dd = dd./max(dd);
    
        if flagMu == 1
            [~,dbest] = max(dd,[],1);
    
        elseif flagMu == 5
            [~,id] = sort(dd,1,'descend');
            dbest  = id(1:min(flagMu,size(dd,1)));
        end
    end

    [frontNo,~] = NDSort([PopObj1(cbest,:);A2Obj],size(A2Obj,1));
    if ~isempty(frontNo(2:end)==2)&&frontNo(1)==1
        PopNew = PopDec1(cbest,:);
    else
        PopNew = PopDec(dbest,:);
    end
    if isempty(PopObj) && isempty(PopObj1)
        PopNew = [];
    end
end