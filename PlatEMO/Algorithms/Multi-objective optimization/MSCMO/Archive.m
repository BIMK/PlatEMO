function Population = Archive(varargin)
% Update the archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select feasible solutions
    [Population,N,index,count] = deal(varargin{1:4});
    if count == 0
        CV = zeros(1,size(Population,2));
    else
        CV = Population.cons;
        CV = CV(:,index(1:count));
    end
    if nargin == 4
        fIndex     = all(CV <= 0,2);
        Population = Population(fIndex);
        if isempty(Population)
            return;
        else
            %% Non-dominated sorting
            [FrontNo,~] = NDSort(Population.objs,1);
            Next = (FrontNo == 1);
            Population = Population(Next);
            Next = true(size(Population,2),1);
            if sum(Next) > N
                Del  = Truncation(Population(Next).objs,sum(Next)-N);
                Temp = find(Next);
                Next(Temp(Del)) = false;
                Population     = Population(Next);      
            end
        end
    else                  
        feasible_number = varargin{5};
        if feasible_number == 0
            return;
        end
        feasible_solutions = Population(1:feasible_number);
        remain_solutions   = Population(feasible_number+1:end);
        [W,~] = UniformPoint(N,size(Population.objs,2));
        itr = 1;
        while size(feasible_solutions,2) < N
            for i = 1 : size(W,1)
                if size(feasible_solutions,2) == N
                    break;
                end
                [~,Region1] = max(1-pdist2(feasible_solutions.objs,W,'cosine'),[],2);
                [~,Region2] = max(1-pdist2(remain_solutions.objs,W,'cosine'),[],2);
                region1     = find(Region1==i);
                region2     = find(Region2==i);
                if length(region1)<itr && ~isempty(region2) 
                    feasible_solutions = [feasible_solutions remain_solutions(region2(1))];
                    remain_solutions(region2(1)) = [];
                end
            end
            itr=itr+1;
        end
        Population=feasible_solutions;
    end  
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
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