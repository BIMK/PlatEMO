function [OffDec,OffMask] = Operator(Problem,ParentDec,ParentMask,rbm,dae,Site,allZero,allOne,FitnessLayer,LayerMax)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [N,D]       = size(ParentDec);
    Parent1Mask = ParentMask(1:N/2,:);
    Parent2Mask = ParentMask(N/2+1:end,:);
    Parent1Dec  = ParentDec(1:N/2,:);
    Parent2Dec  = ParentDec(N/2+1:end,:);
    
    %% Generate mask using crossover (N/2 × D)
    OffMask     = Parent1Mask;
    SparseRate1 = sum(Parent1Mask,2);
    SparseRate2 = sum(Parent2Mask,2);
    Rate        = SparseRate1./(SparseRate1 + SparseRate2);
    
    for i = 1 : N/2
        index = rand(1,D) > Rate(i)/2;
        OffMask(i,index) = Parent2Mask(i,index);
    end
    
    %% Layer-based mutation (still N/2 × D)
    for i = 1 : N/2
        PointUp   = 1;
        PointDown = LayerMax;
        while PointUp < PointDown
            if rand < 0.5
                TargetUpLayer = find(FitnessLayer == PointUp);
                TargetUp      = TargetUpLayer(OffMask(i,TargetUpLayer) == 0);
                if ~isempty(TargetUp) && rand < 0.5
                    TargetUp = datasample(TargetUp,ceil(length(TargetUp)/2));
                    OffMask(i,TargetUp) = 1;
                end
                PointUp = PointUp + 1;
            else
                TargetDownLayer = find(FitnessLayer == PointDown);
                TargetDown      = TargetDownLayer(OffMask(i,TargetDownLayer) == 1);
                if ~isempty(TargetDown) && rand < 0.5
                    TargetDown = datasample(TargetDown,ceil(length(TargetDown)/2));
                    OffMask(i,TargetDown) = 0;
                end
                PointDown = PointDown - 1;
            end
            
            if rand < 0.5 || PointUp >= PointDown
                break;
            end
        end
    end
    
    %% Generate real variables (N/2 × D)
    if any(Problem.encoding~=4)
        if any(Site)
            % Model-based generation
            other = ~allZero & ~allOne;
            if ~isempty(rbm)
                tempMask = OffMask(:,other);
                for i = find(Site)
                    OffTemp = rbm.reduce(tempMask(i,:));
                    OffTemp = rbm.recover(OffTemp);
                    tempMask(i,:) = OffTemp;
                end
                OffMask(:,other)  = tempMask;
                OffMask(:,allOne) = true;
            end
            
            if ~isempty(dae)
                OffDec = zeros(N/2,D);
                for i = find(Site)
                    temp = dae.reduce(Parent1Dec(i,:));
                    temp = dae.recover(temp);
                    OffDec(i,:) = temp;
                end
                % Fill non-Site rows with regular crossover
                nonSite = find(~Site);
                if ~isempty(nonSite)
                    OffDec(nonSite,:) = RealCrossover(Parent1Dec(nonSite,:),Parent2Dec(nonSite,:));
                end
            else
                OffDec = RealCrossover(Parent1Dec,Parent2Dec);
            end
        else
            % Regular crossover
            OffDec = RealCrossover(Parent1Dec,Parent2Dec);
        end
        
        % Mutation and boundary handling
        OffDec = RealMutation(OffDec,Problem.lower,Problem.upper);
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(N/2,D);
    end
end

function Offspring = RealCrossover(Parent1,Parent2)
    [proC,disC] = deal(1,20);
    [N,D]       = size(Parent1);
    
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta = beta.*rand(N,D);
    
    Offspring = (Parent1+Parent2)/2 + beta.*(Parent1-Parent2)/2;
end

function Offspring = RealMutation(Offspring,Lower,Upper)
    [proM,disM] = deal(1,20);
    [N,D] = size(Offspring);
    Lower = repmat(Lower,N,1);
    Upper = repmat(Upper,N,1);
    
    Site = rand(N,D) < proM/D;
    mu   = rand(N,D);
    temp = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    Offspring       = min(max(Offspring,Lower),Upper);
end