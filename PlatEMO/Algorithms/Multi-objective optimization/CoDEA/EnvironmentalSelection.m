function [Population,z,znad] = EnvironmentalSelection(Population,W,N,z,znad,ILid,r)
% The environmental selection of CoDEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    St = find(FrontNo<=MaxFNo);

    %% Normalization
    [PopObj,z,znad] = Normalization(Population(St).objs,z,znad);

    %% CoD-non-dominated sorting
    tFrontNo = CoDSort(PopObj,W,ILid,r);

    %% Selection
    MaxFNo    = find(cumsum(hist(tFrontNo,1:max(tFrontNo)))>=N,1);
    LastFront = find(tFrontNo==MaxFNo);
    LastFront = LastFront(randperm(length(LastFront)));
    tFrontNo(LastFront(1:sum(tFrontNo<=MaxFNo)-N)) = inf;
  
    Next      = St(tFrontNo<=MaxFNo);
    % Population for next generation
    Population = Population(Next);
    

end

function tFrontNo = CoDSort(PopObj,W,ILid,r)
    % Do CoD-non-dominated sorting

    N  = size(PopObj,1);
    [NW,M] = size(W);

    %% Calculate the d1 and d2 values for each solution to each weight
    normP  = sqrt(sum(PopObj.^2,2));
    Cosine = 1 - pdist2(PopObj,W,'cosine');
    d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2);
  
    %% Clustering
    [~,class] = min(d2,[],2);

    %% Sort
    tFrontNo = zeros(1,N);
    
    for i = 1 : NW
        C = find(class==i);
        if(size(C,1) == 0)
            continue;
        end

        if i <= ILid  
            d = max(PopObj(C,:)-W(i,:),[],2);       
            g = ((d+(M)/(1+exp(-1*(M-5.5)*M))*r(i)*((d2(C,i)-min(d2(C,i)))./(max(d2(C,i))-min(d2(C,i))))));
            [~,rank] = sort(g);
        else
             d = 1 - pdist2(PopObj(C,:),ones(1,M)*(1/M),'cosine');
            [~,rank] = sort(d); 
        end
        tFrontNo(C(rank)) = 1 : length(C);
    end



    