function FrontNo = UpdateFront(PopObj,FrontNo,x)
% Update the front No. of each solution when a solution is added or deleted

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(PopObj);
    if nargin < 3
        %% Add a new solution (has been stored in the last of PopObj)
        FrontNo  = [FrontNo,0];
        Move     = false(1,N);
        Move(N)  = true;
        CurrentF = 1;
        % Locate the front No. of the new solution
        while true
            Dominated = false;
            for i = 1 : N-1
                if FrontNo(i) == CurrentF
                    m = 1;
                    while m <= M && PopObj(i,m) <= PopObj(end,m)
                        m = m + 1;
                    end
                    Dominated = m > M;
                    if Dominated
                        break;
                    end
                end
            end
            if ~Dominated
                break;
            else
                CurrentF = CurrentF + 1;
            end
        end
        % Move down the dominated solutions front by front
        while any(Move)
            NextMove = false(1,N);
            for i = 1 : N
                if FrontNo(i) == CurrentF
                    Dominated = false;
                    for j = 1 : N
                        if Move(j)
                            m = 1;
                            while m <= M && PopObj(j,m) <= PopObj(i,m)
                                m = m + 1;
                            end
                            Dominated = m > M;
                            if Dominated
                                break;
                            end
                        end
                    end
                    NextMove(i) = Dominated;
                end
            end
            FrontNo(Move) = CurrentF;
            CurrentF      = CurrentF + 1;
            Move          = NextMove;
        end
    else
        %% Delete the x-th solution
        Move     = false(1,N);
        Move(x)  = true;
        CurrentF = FrontNo(x) + 1;
        while any(Move)
            NextMove = false(1,N);
            for i = 1 : N
                if FrontNo(i) == CurrentF
                    Dominated = false;
                    for j = 1 : N
                        if Move(j)
                            m = 1;
                            while m <= M && PopObj(j,m) <= PopObj(i,m)
                                m = m + 1;
                            end
                            Dominated = m > M;
                            if Dominated
                                break;
                            end
                        end
                    end
                    NextMove(i) = Dominated;
                end
            end
            for i = 1 : N
                if NextMove(i)
                    Dominated = false;
                    for j = 1 : N
                        if FrontNo(j) == CurrentF-1 && ~Move(j)
                            m = 1;
                            while m <= M && PopObj(j,m) <= PopObj(i,m)
                                m = m + 1;
                            end
                            Dominated = m > M;
                            if Dominated
                                break;
                            end
                        end
                    end
                    NextMove(i) = ~Dominated;
                end
            end
            FrontNo(Move) = CurrentF - 2;
            CurrentF      = CurrentF + 1;
            Move          = NextMove;
        end
        FrontNo(x) = [];
    end
end