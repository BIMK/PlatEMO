function Pop = DoubleReproduction(Problem,Pop,GuidingSolution,RefV)
% Generate a promising population by the double reproduction and
% complementary environment selection strategy.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Shufen Qin
% E-mail: shufen.qin@stu.tyust.edu.cn

    %% The First Phase
    OffPopX = GAonce(Problem,Pop.decs,GuidingSolution.decs);
    OffPopX = unique(OffPopX,'rows');
    OffPop  = Problem.Evaluation(OffPopX);
    
    % Environmental Selection
    PopCom = [Pop,GuidingSolution,OffPop];
   
    [associate,Cosinemax] = Assign(PopCom,RefV);
    Ns = length(unique(associate)');
    if Ns < (2/3) * Problem.N
        Pop = DominationSelection(Problem,PopCom);
    else
        Pop = DecompositionSelection(Problem,PopCom,associate,Cosinemax);
    end
    
    %% The Second Phase
    %OffPopX = GA(Pop.decs);
    OffPopX = GATwice(Pop.decs,Problem);
    OffPopX = unique(OffPopX,'rows');
    OffPop  = Problem.Evaluation(OffPopX);
    PopCom  = [Pop,OffPop];
  
    [associate2,Cosinemax2] = Assign(PopCom,RefV);
    Ns = length(unique(associate2));
    if Ns<(2/3)*Problem.N
        Pop = DominationSelection(Problem,PopCom);
    else
        Pop = DecompositionSelection(Problem,PopCom,associate2,Cosinemax2);
    end
end

function [associate,Cosinemax] = Assign(PopCom,RefV)
    % Assign individuals 
    Obj = PopCom.objs;
    Obj = (Obj-repmat(min(Obj),length(PopCom),1))./(repmat(max(Obj),length(PopCom),1)-repmat(min(Obj),length(PopCom),1));
    Cosinetemp = pdist2(Obj,RefV,'cosine');
    Cosine = 1-Cosinetemp;
    [Cosinemax,associate] = max(Cosine,[],2);
end

function Offspring = GAonce(Problem,MatingPool,PopS)
% This function is mainly used in the first reproduction.

    N = size(MatingPool,1);
    RandList = randperm(N);
    MatingPool = MatingPool(RandList, :);

    Ns = size(PopS,1);
    Offspring=[];

    for i = 1: size(MatingPool,1)
        k = randi([1,Ns],1);
        Pop = PopS(k,:);
        Parent = [MatingPool(i,:);Pop];
        Offspringtempt = GAhalfrand(Problem,Parent);
        Offspring = [Offspring;Offspringtempt];
    end
end

function Offspring = GAhalfrand(Problem,Parent)

    %% Parameter setting

	[proC,disC,proM,disM] = deal(0.9,20,1,20);

    Parent1 = Parent(1,:);
    Parent2 = Parent(2,:);
    [N,D]   = size(Parent1);

    %% Genetic operators for real encoding
    % Simulated binary crossover
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    if rand <0.5
        Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;
    else
        Offspring = (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2;
    end
    % Polynomial mutation
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
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

function Offspring = GATwice(MatingPool,Global)
% This function includes the SBX crossover operator and the polynomial
% mutatoion operator.

     MaxOffspring = Global.N;
     [N,D] = size(MatingPool);
     RandList = randperm(N);
     MatingPool = MatingPool(RandList, :);
     if  MaxOffspring < 1 || MaxOffspring > N
         MaxOffspring = N;
     end
     if(mod(N,2) == 1)
         MatingPool = [MatingPool; MatingPool(1,:)];
     end  
     
     ProC = 0.9;
     ProM = 1/D;
     
     DisC = 20;
     DisM = 20;
     Offspring = zeros(N,D);
     %crossover
     for i = 1 : 2 : N
         beta = zeros(1,D);
         miu  = rand(1,D);
         beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(DisC+1));
         beta(miu>0.5)  = (2-2*miu(miu>0.5)).^(-1/(DisC+1));
         beta = beta.*(-1).^randi([0,1],1,D);
         beta(rand(1,D)>ProC) = 1;
         Offspring(i,:)   = (MatingPool(i,:)+MatingPool(i+1,:))/2+beta.*(MatingPool(i,:)-MatingPool(i+1,:))/2;
         Offspring(i+1,:) = (MatingPool(i,:)+MatingPool(i+1,:))/2-beta.*(MatingPool(i,:)-MatingPool(i+1,:))/2;
     end
     Offspring = Offspring(1:MaxOffspring,:);
     
     %mutation
     if MaxOffspring == 1
         MaxValue = Global.upper;
         MinValue = Global.lower;
     else
         MaxValue = repmat(Global.upper,MaxOffspring,1);
         MinValue = repmat(Global.lower,MaxOffspring,1);
     end
     k    = rand(MaxOffspring,D);
     miu  = rand(MaxOffspring,D);
     Temp = k<=ProM & miu<0.5;
     Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Offspring(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
     Temp = k<=ProM & miu>=0.5;
     Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Offspring(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));
     
     Offspring(Offspring>MaxValue) = MaxValue(Offspring>MaxValue);
     Offspring(Offspring<MinValue) = MinValue(Offspring<MinValue);
end