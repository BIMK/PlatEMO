function [OffDec,OffMask] = MSKEA_Operator_sv(Problem,ParentDec,ParentMask,sv)
% The operator of MSKEA guided by pv and sv

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lei Chen
   
    % Parameter setting
    [N,D]       = size(ParentDec);
    if N > 1
        Parent1Mask = ParentMask(1:floor(N/2),:);
        Parent2Mask = ParentMask(floor(N/2)+1:end,:);

        % Crossover for mask
        OffMask = Parent1Mask;
        rate0   = sv;           % The probability that  0 inverts to  1
        rate1   = 1-rate0;      % The probability that  1 inverts to  0

        for i = 1 : floor(N/2)
            diff = find(Parent1Mask(i,:)~=Parent2Mask(i,:));
            temp_rate1=rate1(diff);
            temp_rate0=rate0(diff);
            rate = zeros(1,length(diff));
            rate(logical(OffMask(i,diff)))  = temp_rate1(logical(OffMask(i,diff)));
            rate(logical(~OffMask(i,diff))) = temp_rate0(logical(~OffMask(i,diff)));
            exchange  = rand(1,length(diff)) < rate;
            OffMask(i,diff(exchange))=~OffMask(i,diff(exchange));
        end

        % Mutation for mask
        Mutation_p=1/D;                     % Probability of mutation
        Mu_exchange=rand(floor(N/2),D)<Mutation_p; %The decision variables less than 1/D will be mutated
        for i=1:floor(N/2)
            if sum(Mu_exchange(i,:))
                subscript =find(Mu_exchange(i,:)==1);
                rate=zeros(1,size(subscript,2));
                rate(logical(OffMask(i,subscript)))=rate1(subscript(logical(OffMask(i,subscript))));
                rate(logical(~OffMask(i,subscript)))=rate0(subscript(logical(~OffMask(i,subscript))));
                exchange = rand(1,size(subscript,2)) < rate;
                OffMask(i,subscript(exchange))=~OffMask(i,subscript(exchange));
            end
        end
    else
        % Mutation for mask
        OffMask = ParentMask;
        rate0   = sv;           % The probability that  0 inverts to  1
        rate1   = 1-rate0;      % The probability that  1 inverts to  0
        Mutation_p=1/D;                     % Probability of mutation
        Mu_exchange=rand(1,D)<Mutation_p; %The decision variables less than 1/D will be mutated
        if sum(Mu_exchange)
            subscript =find(Mu_exchange==1);
            rate=zeros(1,size(subscript,2));
            rate(logical(OffMask(:,subscript)))=rate1(subscript(logical(OffMask(:,subscript))));
            rate(logical(~OffMask(:,subscript)))=rate0(subscript(logical(~OffMask(:,subscript))));
            exchange = rand(1,size(subscript,2)) < rate;
            OffMask(:,subscript(exchange))=~OffMask(:,subscript(exchange));
        end
    end
   
    % Crossover and mutation for dec
    if N > 1
        if any(Problem.encoding~=4)
            OffDec = OperatorGAhalf(Problem,ParentDec);
            OffDec(:,Problem.encoding==4) = 1;
        else
            OffDec = ones(floor(N/2),D);
        end
    else
        if any(Problem.encoding~=4)
            OffDec = OperatorMutate(Problem,ParentDec);
            OffDec(:,Problem.encoding==4) = 1;
        else
            OffDec = ParentDec;
        end
    end
end

function Offspring = OperatorMutate(Problem,Parent,Parameter)
%% Mutation for real and integer variables

    if nargin > 2
        [proM,disM] = deal(Parameter{:});
    else
        [proM,disM] = deal(1,20);
    end
    D     = Problem.D;
    Lower = repmat(Problem.lower,1,1);
    Upper = repmat(Problem.upper,1,1);
    Site  = rand(1,D) < proM/D;
    mu    = rand(1,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Parent,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end