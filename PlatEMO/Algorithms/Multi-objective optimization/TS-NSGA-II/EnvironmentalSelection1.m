function Population = EnvironmentalSelection1(OffSpring,W,N,M)
% The selection of TSNSGAII

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    n=2*M;

    ndsPopulation=[];
    Objs=OffSpring.objs;
    [FrontNO,MaxNO]=NDSort(Objs,inf);
    for i=1:MaxNO
        ndsPopulation=cat(2,ndsPopulation,OffSpring(FrontNO==i));
        if length(ndsPopulation)>=N
            break;
        end
    end

    Level=LevelSort(ndsPopulation,n);
    levelPopulation=[];
    for i=1:n
        levelPopulation=cat(2,levelPopulation,ndsPopulation(Level==i));
        %levelPopulation=[levelPopulation;ndsPopulation((Level==i),:)];
        if size(levelPopulation,1)>=N
            break;
        end
    end

    ndsPopulation=levelPopulation;


    Population=[];
    for i=1:size(W,1)
        Zmax=max(ndsPopulation.objs,[],1);
        Zmin=min(ndsPopulation.objs,[],1);
        SPopObj=(ndsPopulation.objs-repmat(Zmin,size(ndsPopulation.objs,1),1))./(repmat(Zmax,size(ndsPopulation.objs,1),1)-repmat(Zmin,size(ndsPopulation.objs,1),1));
        [~,index] = max(1-pdist2(SPopObj,W(i,:),'cosine'));
        Population=[Population,ndsPopulation(index)];
        ndsPopulation=setdiff(ndsPopulation,ndsPopulation(index));
    end
end