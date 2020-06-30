function Population = EnvironmentalSelection(Population,N,TPObj)
% The environmental selection of MaOEA-R&D

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.objs;
    SIn    = find(all(PopObj<=repmat(max(TPObj,[],1),size(PopObj,1),1),2))';
    if length(SIn) > N
        SNd = SIn(NDSort(PopObj(SIn,:),1)==1);
        if length(SNd) > N
            % Select part of the non-dominated solutions
            Next = DiversityOperator(PopObj(SNd,:),N,max([PopObj;TPObj],[],1),min([PopObj;TPObj],[],1));
            Next = SNd(Next);
        else
            % Select several dominated solutions with larger closest
            % distance to non-dominated solutions
            Sd = setdiff(SIn,SNd);
            MinDis = zeros(1,length(Sd));
            for i = 1 : length(Sd)
                MinDis(i) = min(sqrt(sum((repmat(PopObj(Sd(i),:),length(SNd),1)-PopObj(SNd,:)).^2,2)));
            end
            [~,rank] = sort(MinDis,'descend');
            Next = [SNd,Sd(rank(1:N-length(SNd)))];
        end
    else
        % Select several outside solutions which are closest to any target
        % points or their middle points
        TPpair = nchoosek(1:size(TPObj,1),2);
        MiddlePoint = zeros(size(TPpair,1),size(TPObj,2));
        for i = 1 : size(MiddlePoint)
            MiddlePoint(i,:) = (TPObj(TPpair(i,1),:)+TPObj(TPpair(i,2),:))/2;
        end
        Points = [MiddlePoint;TPObj];
        SOut   = setdiff(1:size(PopObj,1),SIn);
        MinDis = zeros(1,length(SOut));
        for i = 1 : length(SOut)
            MinDis(i) = min(sqrt(sum((repmat(PopObj(SOut(i),:),size(Points,1),1)-Points).^2,2)));
        end
        [~,rank] = sort(MinDis);
        Next = [SIn,SOut(rank(1:N-length(SIn)))];
    end
    % Population for next generation
    Population = Population(Next);
end

function Next = DiversityOperator(PopObj,N,fmax,fmin)
% Diversity improvement strategy

    N1   = size(PopObj,1);
    Next = 1 : N1;

    %% Calculate the grid location of each solution
    fmax = repmat(fmax,N1,1);
    fmin = repmat(fmin,N1,1);
    GLoc = N*(PopObj-fmin)./(fmax-fmin);
    
    %% Reduce the number of solutions
    % The grid-based distances between each two solutions
    A = pdist2(GLoc,GLoc);
    A(logical(eye(length(A)))) = inf;
    while length(Next) > N
        [dis,si] = min(A(Next,Next),[],2);
        si = [dis,(1:length(si))',si];
        eliminated = false(1,length(Next));
        while ~isempty(si) && length(Next)-sum(eliminated) > N
            [~,i1] = min(si(:,1));
            i1 = si(i1,2);
            eliminated(i1) = true;
            si(si(:,2)==i1 | si(:,3)==i1,:) = [];
        end
        Next(eliminated) = [];
    end
end