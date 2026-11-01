function Offspring = CompetitiveUpdate(Problem,Loser,Winner,WinnerFitness,operator)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    LoserDec  = Loser.decs;
    WinnerDec = Winner.decs;
    [N,D]     = size(LoserDec);
    LoserVel  = Loser.adds(zeros(N,D));
    WinnerVel = Winner.adds(zeros(N,D));
    r1        = repmat(rand(N,1),1,D);
    r2        = repmat(rand(N,1),1,D);
    OffVel    = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
    OffDec    = LoserDec + OffVel + r1.*(OffVel-LoserVel);
    if operator == 1
        MatingPool      = TournamentSelection(2,N,WinnerFitness);
        WinnerOffspring = OperatorGA(Problem,Winner(MatingPool));
        OffDec          = PolynomialMutation(Problem,OffDec);
        LoserOffspring  = Problem.Evaluation(OffDec,OffVel);
        Offspring       = [LoserOffspring WinnerOffspring];
    else
        [~, minIndex] = min(WinnerFitness);
        WinnerbestDec = WinnerDec(minIndex,:);
        F             = 0.5;
        WinDec        = zeros(N,D);
        for j = 1 : N
            r = randperm(N, 2);
            while any(r == j)
                r = randperm(N, 2);
            end
            WinDec(j,:) = WinnerDec(j,:)+ F.*(WinnerbestDec-WinnerDec(j,:)) + F.*(WinnerDec(r(1),:)-WinnerDec(r(2),:));
        end
        OffDec    = [OffDec;WinDec];
        OffVel    = [OffVel;WinnerVel];
        OffDec    = PolynomialMutation(Problem,OffDec);
        Offspring = Problem.Evaluation(OffDec,OffVel);
    end
end

function OffDec = PolynomialMutation(Problem,OffDec)
    [N,D]  = size(OffDec);
    Lower  = repmat(Problem.lower,N,1);
    Upper  = repmat(Problem.upper,N,1);
    disM   = 20;
    Site   = rand(N,D) < 1/D;
    mu     = rand(N,D);
    temp   = Site & mu<=0.5;
    OffDec = max(min(OffDec,Upper),Lower);
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end