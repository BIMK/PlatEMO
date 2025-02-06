function flag = CalFlag(Population,lastPopulation)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yingwei Li

    FrontNo      = NDSort(Population.objs,inf);
    lastFrontNo  = NDSort(lastPopulation.objs,inf);
    PopObj       = Population.objs;
    lastPopObj   = lastPopulation.objs;
    NDPopObj     = Population(FrontNo==1).objs;
    lastNDPopObj = lastPopulation(lastFrontNo==1).objs;
    lastfmax     = max(lastPopObj(lastFrontNo==1,:),[],1);
    lastfmin     = min(lastPopObj(lastFrontNo==1,:),[],1);
    NormlastND   = (lastNDPopObj-repmat(lastfmin,size(lastNDPopObj,1),1))./repmat(lastfmax-lastfmin,size(lastNDPopObj,1),1);
    lastCon      = sum(NormlastND,2);
    NormND       = (NDPopObj-repmat(lastfmin,size(NDPopObj,1),1))./repmat(lastfmax-lastfmin,size(NDPopObj,1),1);
    Con          = sum(NormND,2);
    HVsum        = CalHV(NDPopObj);
    lastHVsum    = CalHV(lastNDPopObj);
    if lastHVsum < HVsum 
        if min(Con) < min(lastCon)
            flag = 1;
        else
            flag = 2;
        end
    else
        flag = 0;
    end
end

function Score = CalHV(PopObj)
% Calculate the estimated HV value

    RefPoint  = max(PopObj,[],1)*1.1;
    PopObj(any(PopObj>repmat(RefPoint,size(PopObj,1),1),2),:) = [];
    SampleNum = 10000;
    MaxValue  = RefPoint;
    MinValue  = min(PopObj,[],1);
    Samples   = unifrnd(repmat(MinValue,SampleNum,1),repmat(MaxValue,SampleNum,1));
    Domi      = false(1,SampleNum);
    for i = 1 : size(PopObj,1)
        Domi(all(repmat(PopObj(i,:),SampleNum,1)<=Samples,2)) = true;
    end
    Score = prod(MaxValue-MinValue)*sum(Domi)/SampleNum;
end