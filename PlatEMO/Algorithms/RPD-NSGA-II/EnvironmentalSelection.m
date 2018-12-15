function [Population,FrontNo,d2] = EnvironmentalSelection(Population,RPSet,N)
% The environmental selection of RPD-NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Normalization
    PopObj = Population.objs;
    Fmin   = min(PopObj,[],1);
    Fmax   = max(PopObj,[],1);
    PopObj = (PopObj-repmat(Fmin,size(PopObj,1),1))./repmat(Fmax-Fmin,size(PopObj,1),1);
    
    %% Association
    normP   = sqrt(sum(PopObj.^2,2));
    Cosine  = 1 - pdist2(PopObj,RPSet,'cosine');
    d1      = repmat(normP,1,size(RPSet,1)).*Cosine;
    d2      = repmat(normP,1,size(RPSet,1)).*sqrt(1-Cosine.^2);
    [d2,RP] = min(d2,[],2);
    d1      = d1((1:length(RP))'+(RP-1)*length(RP));
    
    %% Favor extreme solutions
    ND              = find(NDSort(PopObj,1)==1);
    [~,Extreme]     = max(PopObj(ND,:),[],1);
    d1(ND(Extreme)) = 0;
    d2(ND(Extreme)) = 0;
    
    %% Non-RPD-dominated sorting
    [FrontNo,MaxFNo] = NRPDDSort(PopObj,d1,d2,RP,N);
    Next = FrontNo < MaxFNo;

    %% Select the solutions in the last front
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(d2(Last));
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    d2         = d2(Next);
end

function [FrontNo,MaxFNo] = NRPDDSort(PopObj,d1,d2,RP,nSort)
% Non-RPD-dominated sorting (based on deductive sort)

    [N,M]   = size(PopObj);
    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while sum(FrontNo<inf) < min(nSort,N)
        MaxFNo = MaxFNo + 1;
        Sorted = FrontNo ~= inf;
        D      = Sorted;
        for i = 1 : N
            if ~D(i)
                for j = i+1 : N
                    if ~D(j)
                        domi = 0;
                        for m = 1 : M
                            if PopObj(i,m) < PopObj(j,m)
                                if domi == -1
                                    domi = 0;
                                    break;
                                else
                                    domi = 1;
                                end
                            elseif PopObj(i,m) > PopObj(j,m)
                                if domi == 1
                                    domi = 0;
                                    break;
                                else
                                    domi = -1;
                                end
                            end
                        end
                        % The definition of dominance is modified here
                        % since the one in the original paper is
                        % questionable
                        if domi == 0 && RP(i)==RP(j)
                            if d1(i)+5*d2(i) < d1(j)+5*d2(j)
                            	domi = 1;
                            elseif d1(i)+5*d2(i) > d1(j)+5*d2(j)
                            	domi = -1;
                            end
                        end
                        if domi == 1
                            D(j) = true;
                        elseif domi == -1
                            D(i) = true;
                            break;
                        end
                    end
                end
                if ~D(i)
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end
end