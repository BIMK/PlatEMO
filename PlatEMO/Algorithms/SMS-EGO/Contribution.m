function fit = Contribution(nleft,PopObj,ypot)
% Calculate the fitness of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    [N,M]      = size(PopObj);
    c          = 1-1/2^M;
    epsilon    = (max(PopObj,[],1)-min(PopObj,[],1))./size(PopObj,1)+c*nleft/500;
    Popepsilon = PopObj-repmat(epsilon,N,1);
    Ny         = size(ypot,1);
    fit        = zeros(Ny,1);
    index      = zeros(Ny,1);
    dindex     = cell(Ny,1);
    refPoint   = max(PopObj,[],1)*1.1;
    HVsum      = CalHV(PopObj,refPoint);
    for i = 1 : Ny
        [index(i),dindex{i}] = Dominance(ypot(i,:),Popepsilon,PopObj);
        if index(i) ==1
            if size(PopObj,2) ==2
                NewObj   = [ypot(i,:);PopObj];
                [~,rank] = sortrows(NewObj);
                j        = find(rank==1);
                if j==1 || j==size(rank,1)
                    domi = any(ypot(i,:)<=refPoint,2) - any(ypot(i,:)>=refPoint,2);
                    if domi == 1
                        fit(i) = -prod(ypot(i,:)-refPoint);
                    else
                        fit(i) = 0;
                    end
                else
                    fit(i) = - (NewObj(rank(j+1),1)-NewObj(rank(j),1)).*(NewObj(rank(j-1),2)-NewObj(rank(j),2));
                end
            else
                fit(i) = HVsum - CalHV([PopObj;ypot(i,:)],refPoint);
            end
        elseif index(i)==2
            domiNo = dindex{i};
            fit(i) = 0;
            for j = 1 : length(domiNo)
                fit(i) = fit(i) -1 + prod(1+(ypot(i,:)-PopObj(domiNo(j),:)));
            end
        else
            domiNo = dindex{i};
            fit(i) = 0;
            flag   = repmat(ypot(i,:),length(domiNo),1)-PopObj(domiNo,:)>0;
            for j = 1 : length(domiNo)
                fit(i) = fit(i) -1 + prod(1+(ypot(i,:)-PopObj(domiNo(j),:)).*flag(j,:));
            end
        end  
    end
end

function [index,dindex] = Dominance(point,Popepsilon,PopObj)
    oobj = repmat(point,size(Popepsilon,1),1);
    domi = any(oobj<=PopObj,2) - any(oobj>=PopObj,2);
    dominated = any(domi==-1);
    if true(dominated)
        index  = 2;
        dindex = find(domi==-1);
    else
        domi = any(oobj<=Popepsilon,2) - any(oobj>=Popepsilon,2);
        dominated = any(domi==-1);
        if true(dominated)
            index  = 3;
            dindex = find(domi==-1);
        else
            index  = 1;
            dindex = 0;
        end
    end
end
 
function Score = CalHV(PopObj,RefPoint)
% Calculate the estimated HV value

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