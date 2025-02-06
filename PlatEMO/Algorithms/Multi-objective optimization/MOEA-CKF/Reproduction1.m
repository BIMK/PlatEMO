function [OffDec,OffMask,len1]=Reproduction1(Problem,Pop1,Dec1,Mask1,FrontNo,CrowdDis,Pop1_Site,Local_Knowlege,Global_Knowlege,NSV,SV,theta,REAL)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    len1 = length(Pop1);
    MatingPool = TournamentSelection(2,len1*2,FrontNo(Pop1_Site),-CrowdDis(Pop1_Site));
    
    % Mask evolution of offspring
    Parent1Mask = Mask1(MatingPool(1:end/2),:);
    Parent2Mask = Mask1(MatingPool((end/2+1):end),:);
    [N,D] = size(Parent1Mask); 
  
    OffMask = false(N,D);
    
    for i = 1 : N
        if rand < theta
            allOne   = Local_Knowlege(1,:);
            other    = Local_Knowlege(3,:);
            TOffMask = false(1,D);          
            TOffMask(other(NSV)) = BinaryCrossover(Parent1Mask(i,other(NSV)),Parent2Mask(i,other(NSV)));
            TOffMask(other(NSV)) = BinaryMutation(TOffMask(other(NSV)));
            
            TOffMask(other(SV)) = BinaryCrossover(Parent1Mask(i,other(SV)),Parent2Mask(i,other(SV)));
            TOffMask(other(SV)) = BinaryMutation(TOffMask(other(SV)));            
            TOffMask(allOne)    = true;  
        else
            allOne   = Global_Knowlege(1,:);
            other    = Global_Knowlege(3,:);
            TOffMask = false(1,D);                                  
            TOffMask(other(NSV)) = BinaryCrossover(Parent1Mask(i,other(NSV)),Parent2Mask(i,other(NSV)));
            TOffMask(other(NSV)) = BinaryMutation(TOffMask(other(NSV)));
            
            TOffMask(other(SV)) = BinaryCrossover(Parent1Mask(i,other(SV)),Parent2Mask(i,other(SV)));
            TOffMask(other(SV)) = BinaryMutation(TOffMask(other(SV)));
            TOffMask(allOne)    = true;  
        end        
        OffMask(i,:) = TOffMask;
    end

    if REAL
        allOne  = Global_Knowlege(1,:);
        allZero = Global_Knowlege(2,:);
        other   = Global_Knowlege(3,:);
        % data preprocessing
        T_Best_Dec       = Dec1(:,~allZero);
        VariableMean     = mean(T_Best_Dec,1);
        VariableNormlize = max(T_Best_Dec)-min(T_Best_Dec)+1e-6;  
        Temp_PopDec      = (T_Best_Dec-repmat(VariableMean,N,1))./repmat(VariableNormlize,N,1);
        % SVD decomposition
        [U,~,~] = svd((1/N).*(Temp_PopDec'*Temp_PopDec));
        % Guide subspace size
        K = sum(other(NSV)) + sum(allOne);
        
        U_reduce      = U(:,1:K);
        Pop_reduce    = Temp_PopDec*U_reduce;
        OffDec_reduce = GAhalfCross(Pop_reduce(MatingPool,:));
        % Restore to original space
        T_OffDec = OffDec_reduce*U_reduce';
        T_OffDec = (T_OffDec.*repmat(VariableNormlize,N,1))+repmat(VariableMean,N,1);
        % Perform polynomial variation
        T_OffDec = PM(Problem,T_OffDec,~allZero);
        %  to reduce the computational complexity of PCA
        if sum(allZero) > 0
            OffDec = zeros(N,D);
            OffDec(:,allZero)  = GAhalf(Problem,Dec1(MatingPool,allZero),allZero);
            OffDec(:,~allZero) = T_OffDec;
        else
            OffDec = T_OffDec;
        end
    else
        OffDec = ones(N,D);
    end
end

function Offspring = PM(Problem,Offspring,site)
    [~,~,proM,disM] = deal(1,20,1,20);
    [N,D] = size(Offspring);
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Lower = Lower(:,site);
    Upper = Upper(:,site);
    
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end

function Offspring = GAhalfCross(Parent)
    [proC,disC,~,~] = deal(1,20,1,20);
    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);
    beta    = zeros(N,D);
    mu      = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;
end

function Offspring = GAhalf(Problem,Parent,allZero)
    [proC,disC,proM,disM] = deal(1,20,1,20);
    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);

    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;

    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Lower = Lower(:,allZero);
    Upper = Upper(:,allZero);

    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end