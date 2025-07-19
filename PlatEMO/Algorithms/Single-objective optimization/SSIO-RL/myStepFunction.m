function [NextObs,Reward,IsDone,LoggedSignals] = myStepFunction(Action,LoggedSignals)
% Custom step function to construct  environment for the function

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Define the environment constants.
    npop = LoggedSignals.npop;
    D    = LoggedSignals.D;
    % Unpack the state vector from the logged signals.
    Lower = repmat(LoggedSignals.Lower,npop,1);
    Upper = repmat(LoggedSignals.Upper,npop,1);
    Pdec  = LoggedSignals.Pdec;
    Pobj  = LoggedSignals.Pobj;

    func    = LoggedSignals.func;
    Reward1 = LoggedSignals.bestObj;
    action  = Action(1:end,:);
    action  = min(max(action,-1),1);

    T  = action(1);
    r  = action(2:21);
    r1 = r(1:10);
    r2 = r(11:end);
    

    WeightP = action(22:end);
    WeightP = min(max(WeightP,1e-10),1);
    WeightProb = WeightP / sum(WeightP);
    Fit     = cumsum(WeightP); 
    Fit     = Fit./max(Fit); % 概率选择分布        

    for k = 1 : LoggedSignals.nOffspring
        Parent1 = Pdec(Tournament(round(max(1,T*LoggedSignals.N)),npop,Pobj),:);
        Parent2 = Pdec(Tournament(1,npop,Pobj),:);
        Parent3 = Pdec(Tournament(1,npop,Pobj),:);
        type    = arrayfun(@(S)find(rand<=Fit,1),1:numel(Parent1));
        type    = reshape(type,size(Parent1));
        Qdec    = Parent1;
        for i = 1 : length(Fit)
            index = type == i;
            Qdec(index) = Parent3(index).*r1(i) + ...
                          Parent2(index).*r2(i) + ...
                          Parent1(index).*(1-r1(i)-r2(i));
        end
        Qdec = min(max(Qdec,Lower),Upper);
       
        Site       = rand(npop,D) < 1/D ;
        mu         = rand(npop,D);
        temp       = Site & mu<=0.5;
        Qdec       = min(max(Qdec,Lower),Upper);
        Qdec(temp) = Qdec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*(1-(Qdec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^21).^(1/21)-1);
        temp       = Site & mu>0.5; 
        Qdec(temp) = Qdec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*(1-(Upper(temp)-Qdec(temp))./(Upper(temp)-Lower(temp))).^21).^(1/21));
        Qdec       = min(max(Qdec,Lower),Upper);
        Qobj       = func(Qdec);
        Pdec       = [Pdec;Qdec];
        Pobj       = [Pobj;Qobj];
        [~,rank]   = sort(Pobj);
        Pdec       = Pdec(rank(1:LoggedSignals.N),:);
        Pobj       = Pobj(rank(1:LoggedSignals.N),:);
    end
    Reward = 10*(Reward1 - Pobj(1))/Reward1;
    if isnan(Reward)
        Reward = 0;
    end
    LoggedSignals.bestObj = Pobj(1);
    LoggedSignals.FES     = LoggedSignals.FES + LoggedSignals.nOffspring * LoggedSignals.npop;
    LoggedSignals.G       = floor(LoggedSignals.FES / LoggedSignals.N);
  
    distance = pdist2(Pdec,Pdec(1,:));
    ccc      = cov(distance,Pobj);
    state1   = ccc(1,2)/(sqrt(ccc(1,1))*sqrt(ccc(2,2)));
    state2   = 1 - (max(Pobj)-min(Pobj))/max(Pobj);
    state3   = 1 - (max(Pobj)-mean(Pobj))/max(Pobj);
    state4   = std(Pobj)/std([repmat(max(Pobj),LoggedSignals.N/2,1);repmat(min(Pobj),LoggedSignals.N/2,1)]);
    state5   = mean(pdist2(Pdec(1,:),Pdec))/pdist2(max(Pdec),min(Pdec));
    state6   =  1 - LoggedSignals.FES/ (LoggedSignals.MaxG*LoggedSignals.N);
    LoggedSignals.Pdec  = Pdec;
    LoggedSignals.Pobj  = Pobj;
    State = [state1;state2;state3;state4;state5;state6];

    LoggedSignals.State = State;
    NextObs = LoggedSignals.State;
    IsDone  = false;

    if LoggedSignals.G>=LoggedSignals.MaxG 
        IsDone = true;
    end

end

function index = Tournament(K,N,varargin)
    if K > 1
        varargin = cellfun(@(S)reshape(S,[],1),varargin,'UniformOutput',false);
        [~,rank] = sortrows([varargin{:}]);
        [~,rank] = sort(rank);
        Parents  = randi(length(varargin{1}),K,N);
        [~,best] = min(rank(Parents),[],1);
        index    = Parents(best+(0:N-1)*K);
    else
        index = randi(length(varargin{1}),1,N);
    end
end
