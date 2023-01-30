function [newMaxP,newMinP,Nonzero] = POSMining(Mask,MaxP,MinP,N)
% Mine the maximum and minimum optimal sparse subspaces

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Nonzero = any(Mask,1);
    maxp    = Mining(Mask(:,Nonzero),MaxP(:,Nonzero),N,true);
    newMaxP = false(size(maxp,1),size(Mask,2));
    newMaxP(:,Nonzero) = maxp;
    newMaxP(1:min(end,size(MaxP,1)),~Nonzero) = MaxP(1:min(end,size(newMaxP,1)),~Nonzero);
    minp    = Mining(Mask(:,Nonzero),MinP(:,Nonzero),N,false);
    newMinP = false(size(minp,1),size(Mask,2));
    newMinP(:,Nonzero) = minp;
    newMinP(1:min(end,size(MinP,1)),~Nonzero) = MinP(1:min(end,size(newMinP,1)),~Nonzero);
    Nonzero = find(Nonzero);
end

function [PopDec,PopObj] = Mining(Mask,PopDec,N,Maximize)
% Evolutionary multi-objective based pattern mining

    %% Population initialization
    Dec = false(N,size(PopDec,2));
    if Maximize
        for i = 1 : N
            Dec(i,:) = any(Mask(randperm(end,ceil(rand.^2*end)),:),1);
        end
    else
        for i = 1 : N
            Dec(i,:) = all(Mask(randperm(end,ceil(rand.^2*end)),:),1);
        end
    end
    PopDec = [PopDec;Dec;xor(Maximize,logical(eye(size(PopDec,2))));Mask];
    Lens   = sum(Mask,2);
    PopObj = CalObj(PopDec,Mask,Lens,Maximize);
	[PopDec,PopObj,FrontNo] = EnvironmentalSelection(PopDec,PopObj,N);
    
    %% Optimization
    for gen = 1 : 10
        MatingPool = TournamentSelection(2,2*N,FrontNo);
        OffDec     = BinaryCrossover(PopDec(MatingPool(1:2:end),:),PopDec(MatingPool(2:2:end),:));
        OffDec     = BinaryMutation(OffDec);
        OffObj     = CalObj(OffDec,Mask,Lens,Maximize);
        [PopDec,PopObj,FrontNo] = EnvironmentalSelection([PopDec;OffDec],[PopObj;OffObj],N);
    end
end

function Objs = CalObj(Decs,T,Lens,Maxmize)
% Calculate the objective values

    [N,D] = size(T);
    Objs  = zeros(size(Decs,1),2);
    Len   = sum(Decs,2);
    if Maxmize
        %% For mining maximum optimal sparse subspaces
        for i = 1 : size(Decs,1)
            Tx = false(1,N);
            for j = 1 : N
                m = 1;
                while m <= D && Decs(i,m) >= T(j,m)
                    m = m + 1;
                end
                Tx(j) = m > D;
            end
            if ~any(Tx)
                Objs(i,:) = 1;
            else
                Objs(i,1) = 1 - mean(Tx);
                Objs(i,2) = 1 - mean(Lens(Tx))./Len(i);
            end
        end
    else
        %% For mining minimum optimal sparse subspaces
        for i = 1 : size(Decs,1)
            Tx = all(T(:,Decs(i,:)),2);
            if ~any(Tx)
                Objs(i,:) = 1;
            else
                Objs(i,1) = 1 - mean(Tx);
                Objs(i,2) = 1 - mean(Len(i)./Lens(Tx));
            end
        end
    end
end