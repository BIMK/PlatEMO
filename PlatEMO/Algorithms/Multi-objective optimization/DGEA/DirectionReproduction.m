function Offspring = DirectionReproduction(Problem,Population,FrontNo,RefNo)
% The direction based offspring generation in DGEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

	[N,D]   = size(Population.decs);
    indexD  = find(FrontNo>1);
    domiN   = numel(indexD);
    indexN  = FrontNo == 1;
	nondec  = Population(indexN).decs;
    domidec = Population(indexD).decs;  
	subN    = max([floor(Problem.N/RefNo),10]);
    Lower   = repmat(Problem.lower,subN,1);
    Upper   = repmat(Problem.upper,subN,1);
    PopDec  = [];
    
    %% Select direction solutions
    startP = nondec(randperm(N-domiN,1),:);
    if domiN < RefNo  
        %This part increases the convergence
        if N <= RefNo
            endP = [domidec;nondec(randperm(N-domiN,N-domiN),:)];
        else
            endP = [domidec;nondec(randperm(N-domiN,RefNo-domiN),:)];
        end
    else
        % This part increase the diversity
        endP = domidec(randperm(domiN,RefNo),:);
    end
    RefNo   = size(endP,1);
    vector  = (endP - repmat(startP,RefNo,1));  
    Direct  = vector./repmat(sum(vector.^2,2).^(1/2),1,D);
    for i = 1 : RefNo
        lambda  = (nondec - repmat(startP,N-domiN,1))*Direct(i,:)';
        sigma   = std(lambda)*(1+((domiN+1)/N)^1);
        OffDec  = repmat(normrnd(0,sigma,[subN,1]),1,D).*repmat(Direct(i,:),subN,1)+repmat(startP,subN,1);
        PopDec  = [PopDec;max(min(OffDec,Upper),Lower)];
    end

    %% Polynomial mutation
	OffDec = PopDec;
    N      = size(OffDec,1);
    Lower  = repmat(Problem.lower,N,1);
    Upper  = repmat(Problem.upper,N,1);
    disM   = 20;
    Site   = rand(N,D) < 1/Problem.D;
    mu     = rand(N,D);
    temp   = Site & mu<=0.5;
    OffDec = max(min(OffDec,Upper),Lower);
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp  = Site & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
	Offspring = Problem.Evaluation(OffDec);
end