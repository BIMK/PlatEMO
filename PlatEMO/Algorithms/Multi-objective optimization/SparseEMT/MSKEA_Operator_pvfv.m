function [OffDec,OffMask] = MSKEA_Operator_pvfv(Problem,ParentDec,ParentMask,pv,fv,delta)
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
        Parent1Mask = ParentMask(1:N/2,:);
        Parent2Mask = ParentMask(N/2+1:end,:);

        % Crossover for mask
        OffMask = Parent1Mask;
        for i = 1 : N/2
            if rand < 0.5
                index = find(Parent1Mask(i,:)&~Parent2Mask(i,:));
                index = index(TS(-pv(index)));
                OffMask(i,index) = 0;
            else
                index = find(~Parent1Mask(i,:)&Parent2Mask(i,:));
                index = index(TS(pv(index)));
                OffMask(i,index) = Parent2Mask(i,index);
            end
        end

        % Mutation for mask
        if rand<(1-delta)
            f_vector=fv;
            f_vector(fv>0)=1;
            for i=1:N/2
                index=find(OffMask(i,:)~=f_vector);
                if rand<0.5
                    index=index(TS(-fv(index)));
                    OffMask(i,index)=1;
                else
                    index=index(TS(fv(index)));
                    OffMask(i,index)=0;
                end
            end
        else
            for i = 1 : N/2
                if rand < 0.5
                    index = find(OffMask(i,:));
                    index = index(TS(-pv(index)));
                    OffMask(i,index) = 0;
                else
                    index = find(~OffMask(i,:));
                    index = index(TS(pv(index)));
                    OffMask(i,index) = 1;
                end
            end
        end
    else
        % Mutation for mask
        OffMask = ParentMask;
        if rand<(1-delta)
            f_vector=fv;
            f_vector(fv>0)=1;
            index=find(OffMask~=f_vector);
            if rand<0.5
                index=index(TS(-fv(index)));
                OffMask(:,index)=1;
            else
                index=index(TS(fv(index)));
                OffMask(:,index)=0;
            end
        else
            if rand < 0.5
                index = find(OffMask);
                index = index(TS(-pv(index)));
                OffMask(:,index) = 0;
            else
                index = find(~OffMask);
                index = index(TS(pv(index)));
                OffMask(:,index) = 1;
            end
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

function index = TS(pv)
% Binary tournament selection

    if isempty(pv)
        index = [];
    else
        index = TournamentSelection(2,1,pv);
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