function [FrontNo,MaxFNo] = NDQSort(varargin)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = varargin{1};
    [N,M]  = size(PopObj);
    if nargin == 2
        nSort = varargin{2};
    else
        PopCon     = varargin{2};
        nSort      = varargin{3};
        Infeasible = any(PopCon>0,2);
        PopObj(Infeasible,:) = repmat(max(PopObj,[],1),sum(Infeasible),1) + repmat(sum(max(0,PopCon(Infeasible,:)),2),1,M);
    end
    if M < 3 || N < 500
        [FrontNo,MaxFNo] = ENS_SS(PopObj,nSort);
    else
        [FrontNo,MaxFNo] = T_ENS(PopObj,nSort);
    end
end

function [FrontNo,MaxFNo] = ENS_SS(PopObj,nSort)
    [PopObj,~,Loc] = unique(PopObj,'rows');
    Table   = hist(Loc,1:max(Loc));
    [N,M]   = size(PopObj);
    Angle   = AssociateAngle(PopObj);
    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
        MaxFNo = MaxFNo + 1;
        for i = 1 : N
            if FrontNo(i) == inf
                Dominated = false;
                for j = i-1 : -1 : 1
                    if FrontNo(j) == MaxFNo
                        m = 2;
                        while m <= M && PopObj(i,m) >= PopObj(j,m)
                            m = m + 1;
                        end
                        if m > M && Angle(i) < Angle(j)
                            Dominated = true;
                        else
                            Dominated = false;
                        end
                        if Dominated || M == 2
                            break;
                        end
                    end
                end
                if ~Dominated
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end

    FrontNo = FrontNo(:,Loc);
end

function [FrontNo,MaxFNo] = T_ENS(PopObj,nSort)
    [PopObj,~,Loc] = unique(PopObj,'rows');
    Table     = hist(Loc,1:max(Loc));
	[N,M]     = size(PopObj);
    FrontNo   = inf(1,N);
    MaxFNo    = 0;
    Forest    = zeros(1,N);
    Children  = zeros(N,M-1);
    LeftChild = zeros(1,N) + M;
    Father    = zeros(1,N);
    Brother   = zeros(1,N) + M;
    [~,ORank] = sort(PopObj(:,2:M),2,'descend');
    ORank     = ORank + 1;
    while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
        MaxFNo = MaxFNo + 1;
        root   = find(FrontNo==inf,1);
        Forest(MaxFNo) = root;
        FrontNo(root)  = MaxFNo;
        for p = 1 : N
            if FrontNo(p) == inf
                Pruning = zeros(1,N);
                q = Forest(MaxFNo);
                while true
                    m = 1;
                    while m < M && PopObj(p,ORank(q,m)) >= PopObj(q,ORank(q,m))
                        m = m + 1;
                    end
                    if m == M
                        break;
                    else
                        Pruning(q) = m;
                        if LeftChild(q) <= Pruning(q)
                            q = Children(q,LeftChild(q));
                        else
                            while Father(q) && Brother(q) > Pruning(Father(q))
                                q = Father(q);
                            end
                            if Father(q)
                                q = Children(Father(q),Brother(q));
                            else
                                break;
                            end
                        end
                    end
                end
                if m < M
                    FrontNo(p) = MaxFNo;
                    q = Forest(MaxFNo);
                    while Children(q,Pruning(q))
                        q = Children(q,Pruning(q));
                    end
                    Children(q,Pruning(q)) = p;
                    Father(p) = q;
                    if LeftChild(q) > Pruning(q)
                        Brother(p)   = LeftChild(q);
                        LeftChild(q) = Pruning(q);
                    else
                        bro = Children(q,LeftChild(q));
                        while Brother(bro) < Pruning(q)
                            bro = Children(q,Brother(bro));
                        end
                        Brother(p)   = Brother(bro);
                        Brother(bro) = Pruning(q);
                    end
                end
            end
        end
    end
    FrontNo = FrontNo(:,Loc);
end