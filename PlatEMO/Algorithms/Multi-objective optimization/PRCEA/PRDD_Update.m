function [return_pop,return_Fitness] = PRDD_Update(mainPop,Population,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    return_pop     = [];
    return_Fitness = [];
    maxfit         = 0;
    input_cons     = mainPop.cons;
    input_cons(input_cons<0) = 0;
    input_cons = sum(input_cons,2);

    S1 = [];
    S2 = [];
    S3 = [];
    if sum(input_cons<=0) == 0
        S1 = [S1,Population];
    else
        FP = mainPop(input_cons<=0);
        [FrontNo,~] = NDSort(FP.objs,1);
        FP = FP(FrontNo ==1);
        for i = 1 : length(Population)
            Diff = sign(Population(i).objs-FP.objs);
            if sum(max(Diff,[],2)-min(Diff,[],2)==2)==length(FP)
                S1 =[S1,Population(i)];
            elseif any(sum(Diff<=0,2)==size(FP.objs,2))
                S2 = [S2,Population(i)];
            elseif any(sum(Diff>=0,2)==size(FP.objs,2))
                S3 = [S3,Population(i)];
            end
        end
    end

    if length(S1)>=N
        %EPD for S1
        [return_pop,return_Fitness] = CDP_EPD_Update(S1,N,false);
    else
        if ~isempty(S1)
            return_pop     = S1;
            return_Fitness = CalFitness([S1.objs,sum(max(S1.cons,0),2)]);
            maxfit         = max(return_Fitness);
        end
        %ERPD for S2
        if length(S2) >= N-length(S1)
            Fitness = CalFitness([-S2.objs,sum(max(S2.cons,0),2)]);
            Next    = Fitness < 1;
            if sum(Next) <= N - length(S1)
                [~,Rank] = sort(Fitness);
                Next(Rank(1:N - length(S1) )) = true;
            elseif sum(Next) > N - length(S1)
                Del  = Truncation(S2(Next).objs, sum(Next)-(N-length(S1)));
                Temp = find(Next);
                Next(Temp(Del)) = false;
            end
            S2             = S2(Next);
            Fitness        = Fitness(Next) + maxfit;      
            return_pop     = [return_pop,S2];
            return_Fitness = [return_Fitness,Fitness];
        %PD for S3
        else
            if ~isempty(S2)
                return_pop     = [return_pop,S2];
                Fitness        = CalFitness([-S2.objs,sum(max(S2.cons,0),2)])+maxfit;
                return_Fitness = [return_Fitness,Fitness];
                maxfit         = max(return_Fitness);
            end

            Fitness = CalFitness(S3.objs); 
            Next    = Fitness < 1;
            if sum(Next) <= N - length(return_pop)
                [~,Rank] = sort(Fitness);
                Next(Rank(1:N - length(return_pop) )) = true;
            elseif sum(Next) > N - length(return_pop)
                Del  = Truncation(S3(Next).objs, sum(Next)-(N-length(return_pop)));
                Temp = find(Next);
                Next(Temp(Del)) = false;
            end
            S3             = S3(Next);
            Fitness        = Fitness(Next) + maxfit;      
            return_pop     = [return_pop,S3];
            return_Fitness = [return_Fitness,Fitness];
        end
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