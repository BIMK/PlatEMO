function Offspring01 = Adaptive_sampling(CAobj,DAobj,CAdec,DAdec,DAvar,DA,CA1,mu,p,phi)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhenshou Song

    Ideal_Point = min([CAobj;DAobj],[],1);

    if  size(DAdec,1) <= mu
        flag = 1;
    else
        flag = Cal_Convergence(CAobj,DAobj,Ideal_Point);
    end

    if flag == 1
        % convergence sampling strategy
        N = size(CAobj,1);
        CAobj01 = (CAobj-repmat(min(CAobj),N,1))./(repmat(max(CAobj)-min(CAobj),N,1));
        I = zeros(N);
        for i = 1:N
            for j = 1:N
                I(i,j) = max(CAobj01(i,:)-CAobj01(j,:));
            end
        end
        C = max(abs(I));
        F = sum(-exp(-I./repmat(C,N,1)/0.05)) + 1;
        Choose = 1:N;
        while length(Choose) > mu
            [~,x] = min(F(Choose));
            F = F + exp(-I(Choose(x),:)/C(Choose(x))/0.05);
            Choose(x) = [];
        end
        Offspring01 = CAdec(Choose,:);
    else
        if size(DAdec,1) <= mu
            Offspring01 = DAdec;
        else
            if PD(DAobj) < PD(DA.objs)
                % uncertainty sampling strategy
                An = size(DAvar,1);
                Choose = zeros(1,5);
                for i = 1:mu
                    A_num = randperm(size(DAvar,1));
                    Uncertainty = mean(DAvar(A_num(1:ceil(phi*An)),1:size(DAobj,2)),2);
                    [~,best]    = max(Uncertainty);
                    Choose (i)     = A_num(best);
                end
                Offspring01 = DAdec(Choose ,:);
            else
                % diversity sampling strategy
                DA_Nor = (DA.objs - repmat(min([DAobj;DA.objs],[],1),length(DA),1))...
                    ./repmat(max([DAobj;DA.objs],[],1) - min([DAobj;DA.objs],[],1),length(DA),1);
                DA_Nor_pre = (DAobj - repmat(min([DAobj;DA.objs],[],1),size(DAobj,1),1))...
                    ./repmat(max([DAobj;DA.objs],[],1) - min([DAobj;DA.objs],[],1),size(DAobj,1),1);
                N  = size(DA_Nor,1);
                Pop = [DA_Nor;DA_Nor_pre];
                Pop_dec = [DA.decs;DAdec];
                NN = size(Pop,1);
                Choose = false(1,NN);
                Choose(1:N) = true;
                MaxSize = N+mu;
                Distance = inf(N);
                for i = 1 : NN-1
                    for j = i+1 : NN
                        Distance(i,j) = norm(Pop(i,:)-Pop(j,:),p);
                        Distance(j,i) = Distance(i,j);
                    end
                end
                Offspring01 = [];
                while sum(Choose) < MaxSize
                    Remain = find(~Choose);
                    [~,x]  = max(min(Distance(~Choose,Choose),[],2));
                    Choose(Remain(x)) = true;
                    Offspring01 = [Offspring01;Pop_dec(Remain(x),:)];
                end
            end
        end
    end
end

function Score = PD(PopObj)
% Pure diversity

    N = size(PopObj,1);
    C = false(N);
    C(logical(eye(size(C)))) = true;
    D = pdist2(PopObj,PopObj,'minkowski',0.1);
    D(logical(eye(size(D)))) = inf;
    Score = 0;
    for k = 1 : N-1
        while true
            [d,J] = min(D,[],2);
            [~,i] = max(d);
            if D(J(i),i) ~= -inf
                D(J(i),i) = inf;
            end
            if D(i,J(i)) ~= -inf
                D(i,J(i)) = inf;
            end
            P = any(C(i,:),1);
            while ~P(J(i))
                newP = any(C(P,:),1);
                if P == newP
                    break;
                else
                    P = newP;
                end
            end
            if ~P(J(i))
                break;
            end
        end
        C(i,J(i)) = true;
        C(J(i),i) = true;
        D(i,:)    = -inf;
        Score     = Score + d(i);
    end
end