function Population = EnvironmentalSelection(Population,K,pc)
% The environmental selection of SRA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N      = length(Population);
    PopObj = Population.objs;
    
    %% Compute indicator values of I1 (epsilon+)
    I = zeros(N);
    for i = 1 : N
        for j = 1 : N
            I(i,j) = max(PopObj(i,:)-PopObj(j,:));
        end
    end
    I1 = sum(-exp(-I./0.05)) + 1;
    
    %% Compute indicator values of I2 (SDE)
    Distance = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = 1 : i-1
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    I2 = min(Distance,[],2);
    
    %% Stochastic ranking based selection
    Rank = 1 : N;
    for sweepCounter = 1 : ceil(N/2)
        swapdone = false;
        for j = 1 : N-1
            if rand < pc
                if I1(Rank(j)) < I1(Rank(j+1))
                    Rank     = Rank([1:j-1,j+1,j,j+2:end]);
                    swapdone = true;
                end
            else
                if I2(Rank(j)) < I2(Rank(j+1))
                    Rank     = Rank([1:j-1,j+1,j,j+2:end]);
                    swapdone = true;
                end
            end
        end
        if ~swapdone
            break;
        end
    end
    Population = Population(Rank(1:K));
end