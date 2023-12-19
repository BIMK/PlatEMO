function [Input,Output,Pa,Pmid] = CalFitnessPC(PopObj,PopDec,rate)
% Calculate the fitness value by using a new way  considering 
% convergence and the diversity 

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N      = size(PopObj,1);
    Zmin   = min(PopObj,[],1);
    Zmax   = max(PopObj);            
    PopObj = (PopObj-repmat(Zmin,N,1))./(repmat(Zmax,N,1)-repmat(Zmin,N,1));
    SDE    = zeros(N,1);
    for i = 1 : N
        SPopuObj = PopObj;
        Temp     = repmat(PopObj(i,:),N,1);
        Shifted  = PopObj < Temp;
        SPopuObj(Shifted) = Temp(Shifted);                                    
        Distance  = pdist2(real(PopObj(i,:)),real(SPopuObj));
        [~,index] = sort(Distance,2);
        Dk = Distance(index(floor(sqrt(N))+1)); % Dk denotes the distance of solution i and its floor(sqrt(N)+1)-th nearest neighbour
        SDE(i)=2./(Dk+2);
    end
    Objs = SDE;
    
   %% Detect the dominance relation between each two solutions
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N                 
            k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end

   %% Calculate S(i)
    S = sum(Dominate,2);

    %% Calculate R(i)
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end

   %% Calculate D(i)
    Distance = pdist2(real(Objs),real(Objs),'cosine');
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D = 1./(Distance(:,floor(sqrt(N)))+2);

   %% Calculate the fitnesses
    Rmin = min(R);
    Rmax = max(R);
    R = (R-repmat(Rmin,1,N))./(repmat(Rmax,1,N)-repmat(Rmin,1,N));
    R = R';
    Fitness =  rate*R + (1-rate)*D;

    [~, index] = sort(Fitness);
    Input  = [PopDec(index(1:ceil(N/4)),:);PopDec(index(end-(ceil(N/2)-ceil(N/4))+1:end),:)];
    Input  = [Input,[Fitness(index(1:ceil(N/4)));Fitness(index(end-(ceil(N/2)-ceil(N/4))+1:end))]];
    Output = zeros(ceil(N/2),1);
    Output(1:ceil(N/4))     = 2;
    Output(ceil(N/4)+1:end) = 1;
    Pa   = Input(1:ceil(N/4),:);
    Pmid = Input(ceil(N/4)-floor(N/8):ceil(N/4)+floor(N/8),:);
end