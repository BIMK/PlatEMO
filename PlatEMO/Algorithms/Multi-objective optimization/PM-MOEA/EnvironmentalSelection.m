function varargout = EnvironmentalSelection(varargin)
% The environmental selection of PM-MOEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete duplicated solutions and non-dominated sorting
    if nargin == 3
        [PopDec,PopObj,N] = deal(varargin{:});
        [PopDec,uni] = unique(PopDec,'rows','stable');
        PopObj = PopObj(uni,:);
        [FrontNo,MaxFNo] = NDSort(PopObj,N);
    elseif nargin == 4
        [Population,Dec,Mask,N] = deal(varargin{:});
        PopObj = Population.objs;
        [FrontNo,MaxFNo] = NDSort(PopObj,Population.cons,N);
    end
    Next = FrontNo <= MaxFNo;
    
    %% Truncate the solutions in the last front
    Last = find(FrontNo==MaxFNo);
    if nargin == 3
        Del = Truncation(double(PopDec(Last,:)),sum(Next)-N);
    elseif nargin == 4
        Del = Truncation(PopObj(Last,:),sum(Next)-N);
    end
    Next(Last(Del)) = false;
    
    %% Population for next generation
    FrontNo = FrontNo(Next);
    if nargin == 3
        PopDec    = PopDec(Next,:);
        PopObj    = PopObj(Next,:);
        varargout = {PopDec,PopObj,FrontNo};
    elseif nargin == 4
        Population = Population(Next);
        Dec        = Dec(Next,:);
        Mask       = Mask(Next,:);
        varargout  = {Population,Dec,Mask,FrontNo};
    end
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end