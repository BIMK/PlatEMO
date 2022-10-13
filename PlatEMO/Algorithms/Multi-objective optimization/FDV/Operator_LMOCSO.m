function [OffDec,OffVel] = Operator_LMOCSO(Problem,Loser,Winner,Rate)
% The competitive swarm optimizer of LMOCSO

%  Copyright (C) 2021 Xu Yang
%  Xu Yang <xuyang.busyxu@qq.com> or <xuyang369369@gmail.com>

    %% Parameter setting
    LoserDec  = Loser.decs;
    WinnerDec = Winner.decs;
    [N,D]     = size(LoserDec);
	LoserVel  = Loser.adds(zeros(N,D));
    WinnerVel = Winner.adds(zeros(N,D));
    
    %% Competitive swarm optimizer
    r1     = repmat(rand(N,1),1,D);
    r2     = repmat(rand(N,1),1,D);
    
    OffVel = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
    OffDec = LoserDec + OffVel;
    
    if Problem.FE/Problem.maxFE < Rate
        LoserVel1 = rand(N,D);
        OffVel1 = r1.*LoserVel1 + r2.*(WinnerDec-LoserDec);
        OffDec1 = LoserDec + OffVel1 + r1.*(OffVel1-LoserVel1);
        
        OffDec = [OffDec;OffDec1];
        OffVel = [OffVel;OffVel1];
    end
    
    %% Add the winners
    OffDec = [OffDec;WinnerDec];
    OffVel = [OffVel;WinnerVel];
 
    %% Polynomial mutation
    [N,D] = size(OffDec);
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    disM  = 20;
    Site  = rand(N,D) < 1/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    OffDec       = max(min(OffDec,Upper),Lower);
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp  = Site & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));          
end