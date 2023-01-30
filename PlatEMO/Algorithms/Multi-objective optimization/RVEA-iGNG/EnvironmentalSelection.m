function [Population,net,V,Archive,scale,genFlag] = EnvironmentalSelection(Population,V,theta,net,params,Archive,Problem,scale,zmin,genFlag)
% The environmental selection of RVEA

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu

    Population = Population(NDSort(Population.objs,1)==1);
    PopObj = Population.objs;
    [N,M]  = size(PopObj);
    %% Translate the population
    PopObj = PopObj - repmat(zmin,N,1);
    %% delete the outliers
    d = sqrt(sum(PopObj.^2,2));
    meanD = sum(d,1)/size(PopObj,1);
    delete = find(d>10*meanD);
    Population(delete)=[];
    PopObj(delete,:)=[];
    SavePopObj = Population.objs;


    Archive = UpdateArchive(Population,Archive,2*Problem.N);
    ArcObj = Archive.objs;

    wholeObj = [SavePopObj;ArcObj];
    Population = [Population Archive];
    [c,ia,ic] = unique(wholeObj,'rows');
    Population = Population(ia);
    wholeObj = wholeObj(ia,:);
    wholeObj = wholeObj - repmat(zmin,size(wholeObj,1),1);



    wholeObj1 = wholeObj./scale;
    temp1 = wholeObj1./sum(wholeObj1,2);


    PopObj = wholeObj;
    [N,M]  = size(PopObj);
    fr = 0.1;

    gen    = ceil(Problem.FE/Problem.N);
    maxgen = ceil(Problem.maxFE/Problem.N);
    if ~mod(gen,ceil(fr*maxgen))&& gen <= round(1*maxgen)
        scale = max(ArcObj,[],1)-min(ArcObj,[],1);
    end

    if size(temp1,1) > 2&&isempty(find(isnan(temp1)==true))&& gen <= round(1*maxgen) && isempty(genFlag)
        [V,net,genFlag] = TrainGrowingGasNet(V,temp1,net,scale,params,Problem,[[];ArcObj],genFlag,zmin);
    end
    NV     = size(V,1);

    %% Calculate the degree of violation of each solution

    CV = sum(max(0,Population.cons),2);

    %% Calculate the smallest angle value between each vector and others
    cosine = 1 - pdist2(V,V,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2);

    %% Associate each solution to a reference vector
    Angle = acos(1-pdist2(PopObj,V,'cosine'));
    [~,associate] = min(Angle,[],2);

    %% Select one solution for each reference vector
    Next = zeros(1,NV);
    for i = unique(associate)'
        current1 = find(associate==i & CV==0);
        current2 = find(associate==i & CV~=0);
        if ~isempty(current1)
            % Calculate the APD value of each solution
            APD = (1+M*theta*Angle(current1,i)/gamma(i)).*sqrt(sum(PopObj(current1,:).^2,2));
            % Select the one with the minimum APD value
            [~,best] = min(APD);
            Next(i)  = current1(best);
        elseif ~isempty(current2)
            % Select the one with the minimum CV value
            [~,best] = min(CV(current2));
            Next(i)  = current2(best);
        end
    end
    % Population for next generation
    Population1 = Population(Next(Next~=0));
    %% select the corner solutions in each generation

    fm = [];
    selectedFirst = unique([Next(Next~=0) fm]);
    Population = [Population(selectedFirst)];

    if length(Population) > Problem.N 
        PopObj = Population.objs;
        zmax = max(PopObj,[],1);
        PopObj = PopObj - repmat(zmin,size(PopObj,1),1);

        temp1 = PopObj;

        Choose = false(1,size(temp1,1));
        [~,Extreme1] = min(temp1,[],1);
        [~,Extreme2] = max(temp1,[],1);
        Choose(Extreme1) = true;
        Choose(Extreme2) = true;  

        while sum(Choose) < Problem.N
            ind = find(Choose== false);
            choId = setdiff(1:size(temp1,1),ind); 
            PopObj1temp = temp1(choId,:);     
            WholeObjtemp = temp1(ind,:);
            dis =  pdist2(WholeObjtemp,PopObj1temp);
            [mindis,] = min(dis,[],2);
            [~,associate] = max(mindis,[],1);
            Choose(ind(associate)) = true;
        end
        Population = Population(Choose);
    end
end