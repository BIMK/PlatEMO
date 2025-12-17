classdef Block_Selection < BLOCK
% Environmental selection
% nSolutions --- 100 --- Number of retained solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        nSolutions;     % <hyperparameter> Number of retained solutions
    end
    methods
        %% Default settings of the block
        function obj = Block_Selection(nSolutions)
            obj.nSolutions = nSolutions;    % Number of retained solutions
        end
        %% Main procedure of the block
        function Main(obj,Problem,Precursors,Ratio)
            Population = obj.Gather(Problem,Precursors,Ratio,1,1);
            if Problem.M == 1	% For single-objective optimization
                [~,rank]   = sort(FitnessSingle(Population));
                obj.output = Population(rank(1:min(end,obj.nSolutions)));
            else                % For multi- and many-objective optimization
                [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,obj.nSolutions);
                Next = find(FrontNo<MaxFNo);
                Last = find(FrontNo==MaxFNo);
                Del  = Truncation(Population(Last).objs,Population([Last,Next]).objs,length([Last,Next])-obj.nSolutions);
                obj.output = Population([Next,Last(~Del)]);
            end
        end
    end
end

function Del = Truncation(PopObjLast,PopObjAll,K)
    if size(PopObjLast,2) < 4
        Distance = pdist2(PopObjLast,PopObjAll);
        Distance(logical(eye(size(PopObjLast,1)))) = inf;
    else
        Distance = inf(size(PopObjLast,1),size(PopObjAll,1));
        for i = 1 : size(PopObjLast,1)
            SPopObj = max(PopObjAll,repmat(PopObjLast(i,:),size(PopObjAll,1),1));
            for j = [1:i-1,i+1:size(PopObjAll,1)]
                Distance(i,j) = norm(PopObjLast(i,:)-SPopObj(j,:));
            end
        end
    end
    Del = false(1,size(PopObjLast,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,[Remain,size(PopObjLast,1)+1:end]),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end